"""
ARIADNE - Plots module

Make figures.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import terms





def dssr_analysis(dssr_info):
    for section, section_info in dssr_info.items():
        print(section)
        sns.set_style("darkgrid")
        sns.despine()
        # bond lengths vs bond type
        fig = plt.figure()
        design_order = ["LibFig_rPB66", "LibFig_rOct66", "rPB66_v2", "rOct66_v2"]
        ax = sns.stripplot(data=section_info, x="parameter", y="value", hue="design", hue_order=design_order, orient="v", size=2, dodge=True)
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        aspect = (x_max - x_min) / (y_max - y_min)
        ax.set_aspect(aspect)
        fig.set_size_inches(10, 6)
        #ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
        #ax.grid(axis="y", which="minor")
        plt.legend(bbox_to_anchor=(1.0, 1.0))
        plt.title(section)
        plt.savefig(f"dssr_{section.replace(' ', '_')}.png", dpi=200)
        plt.close()


def bond_length_distribution(bond_lengths, design_order=None):
    sns.set_style("darkgrid")
    sns.despine()
    # bond lengths vs bond type
    fig = plt.figure()
    fig.set_size_inches(14, 4)
    ax = sns.stripplot(data=bond_lengths, x="bond length (Å)", y="bond type", hue="design", hue_order=design_order, orient="h", size=2, dodge=True)
    ax.set_aspect(1.0)
    ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
    ax.grid(axis="y", which="minor")
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.title("Length of bond vs. location of bond across origami designs")
    plt.savefig("bond_length_vs_type.png", dpi=200)
    plt.close()
    # length of O3'-P bond vs location of base
    fig = plt.figure()
    fig.set_size_inches(14, 8)
    bond_lengths_inter = bond_lengths.loc[bond_lengths["bond type"] == "O3'-P length"]
    order = sorted(set(bond_lengths_inter["Feature"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="Feature", hue="design", order=order, hue_order=design_order, orient="h", size=2, dodge=True)
    ax.set_aspect(0.5)
    ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
    ax.grid(axis="y", which="minor")
    plt.legend(bbox_to_anchor=(1.01, 1.0))
    plt.title("Length of O3'-P bond vs. location of base across origami designs")
    plt.savefig("bond_length_vs_location.png", dpi=200)
    plt.close()
    # axial length of O3'-P bond vs location of base
    fig = plt.figure()
    fig.set_size_inches(14, 8)
    bond_lengths_inter = bond_lengths.loc[bond_lengths["bond type"] == "O3'-P axial"]
    order = sorted(set(bond_lengths_inter["Feature"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="Feature", hue="design", order=order, hue_order=design_order, orient="h", size=2, dodge=True)
    ax.set_aspect(0.5)
    ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
    ax.grid(axis="y", which="minor")
    plt.legend(bbox_to_anchor=(1.01, 1.0))
    plt.title("Axial length of O3'-P bond vs. location of base across origami designs")
    plt.savefig("bond_length_axial_vs_location.png", dpi=200)
    plt.close()
    # planar length of O3'-P bond vs location of base
    fig = plt.figure()
    fig.set_size_inches(14, 8)
    bond_lengths_inter = bond_lengths.loc[bond_lengths["bond type"] == "O3'-P planar"]
    order = sorted(set(bond_lengths_inter["Feature"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="Feature", hue="design", order=order,
                       hue_order=design_order, orient="h", size=2, dodge=True)
    ax.set_aspect(0.5)
    ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
    ax.grid(axis="y", which="minor")
    plt.legend(bbox_to_anchor=(1.01, 1.0))
    plt.title("Planar length of O3'-P bond vs. location of base across origami designs")
    plt.savefig("bond_length_planar_vs_location.png", dpi=200)
    plt.close()


def cmap_tricolor(value):
    COLORS = ["#7fc2ff", 0.01, "#d1bb3a", 0.03, "#c61e61"]
    if np.isnan(value):
        return "#e0e0e0"
    assert 0 <= value <= 1
    color = COLORS[0]
    for item in COLORS[1:]:
        if isinstance(item, float):
            if value <= item:
                return color
        else:
            color = item
    return color


def cmap_named(name, value):
    if np.isnan(value):
        return "#e0e0e0"
    assert 0 <= value <= 1
    color = plt.get_cmap(name)(value)
    return color


def secondary_structure_signal(fname, edges, g_up, g_dn, g_ax, base_info, signals):
    """
    Plot the secondary structure of an origami, optionally overlaid with a chemical probing signal.
    :return:
    """
    X_INC = 1  # horizontal space between adjacent bases
    Y_INC = 12  # vertical space between adjacent edges
    Y_SEP = 4  # vertical space between helices within an edge
    NUM_SEP = 3  # vertical space between helix and base number
    NUM_PERIOD = 20  # frequency of numbered bases
    BASE_SIZE = 50  # area of the base
    TM_HEIGHT = 3  # height of terminus markers
    SEQ_PARAMS = {"ha": "center", "va": "center", "size": 6}
    NUM_PARAMS = {"ha": "center", "va": "center", "size": 8}
    SCAF_XO_COLOR = "#000000"
    STAP_XO_COLOR = "#50e350"
    SCAF_XO = f"{terms.SCAF}_{terms.SCAF_XO}"
    STAP_XO = f"{terms.SCAF}_{terms.STAP_XO}"
    SCAF_TM = f"{terms.SCAF}_{terms.SCAF_TM}"
    STAP_TM = f"{terms.SCAF}_{terms.STAP_TM}"
    zorder = {"feature": 1, "Base": 2, "text": 3}
    seqs = dict(zip(base_info["CanDo number"], base_info["Base"]))
    locs = dict(zip(base_info["CanDo number"], base_info["Feature"]))
    xs = dict()
    ys = dict()
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    plt.axis("off")  # hide axes and ticks
    for edge_i, double_helices in enumerate(edges):
        for helix_i, helix in enumerate(double_helices):
            for pos, num in enumerate(helix):
                x = pos * X_INC
                y = -edge_i * Y_INC + Y_SEP * [0, -1][helix_i]
                xs[num], ys[num] = x, y
                plt.scatter([x], [y], c=[cmap_named("inferno", signals.get(num, np.nan))], s=BASE_SIZE, zorder=zorder["Base"])
                plt.text(x, y, seqs[num], zorder=zorder["text"], c="#ffffff", **SEQ_PARAMS)
                if num == 1 or num % NUM_PERIOD == 0:
                    plt.text(x, y + NUM_SEP * [1, -1][helix_i], num, **NUM_PARAMS)
                loc = locs[num]
                if loc.startswith(SCAF_XO):
                    # scaffold crossovers
                    if loc.endswith("5"):
                        # x position of the base to which it is crossing over
                        x_xo = xs.get(g_dn[num])
                        y_xo = ys.get(g_dn[num])
                    elif loc.endswith("3"):
                        x_xo = xs.get(g_up[num])
                        y_xo = ys.get(g_up[num])
                    else:
                        raise ValueError()
                    if x_xo is not None:
                        plt.plot([x, x_xo], [y, y_xo], c=SCAF_XO_COLOR, zorder=zorder["feature"])
                elif loc.startswith(STAP_XO):
                    # staple crossovers
                    if loc.endswith("5"):
                        # x position of the base to which it is crossing over
                        x_xo = xs.get(g_ax[g_up[g_ax[num]]])
                        y_xo = ys.get(g_ax[g_up[g_ax[num]]])
                    elif loc.endswith("3"):
                        # x position of the base to which it is crossing over
                        x_xo = xs.get(g_ax[g_dn[g_ax[num]]])
                        y_xo = ys.get(g_ax[g_dn[g_ax[num]]])
                    else:
                        raise ValueError()
                    if x_xo is not None:
                        plt.plot([x, x_xo], [y, y_xo], c=STAP_XO_COLOR, zorder=zorder["feature"])
                elif loc.startswith(SCAF_TM):
                    # staple crossovers
                    if loc.endswith("5"):
                        x_tm = x + X_INC / 2 * [-1, 1][helix_i]
                    elif loc.endswith("3"):
                        x_tm = x + X_INC / 2 * [1, -1][helix_i]
                    else:
                        raise ValueError()
                    if x_tm is not None:
                        plt.plot([x_tm, x_tm], [y, y + TM_HEIGHT * [1, -1][helix_i]], c=SCAF_XO_COLOR, zorder=zorder["feature"])
                elif loc.startswith(STAP_TM):
                    # staple crossovers
                    if loc.endswith("5"):
                        x_tm = x + X_INC / 2 * [1, -1][helix_i]
                    elif loc.endswith("3"):
                        x_tm = x + X_INC / 2 * [-1, 1][helix_i]
                    else:
                        raise ValueError()
                    if x_tm is not None:
                        plt.plot([x_tm, x_tm], [y, y + TM_HEIGHT * [1, -1][helix_i]], c=STAP_XO_COLOR, zorder=zorder["feature"])
    plt.savefig(fname)
