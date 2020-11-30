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


def bond_length_distribution(bond_lengths):
    sns.set_style("darkgrid")
    sns.despine()
    # bond lengths vs bond type
    fig = plt.figure()
    fig.set_size_inches(14, 4)
    design_order = ["LibFig_rPB66", "LibFig_rOct66", "rPB66_v2", "rOct66_v2"]
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
    order = sorted(set(bond_lengths_inter["location"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="location", hue="design", order=order, hue_order=design_order, orient="h", size=2, dodge=True)
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
    order = sorted(set(bond_lengths_inter["location"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="location", hue="design", order=order, hue_order=design_order, orient="h", size=2, dodge=True)
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
    order = sorted(set(bond_lengths_inter["location"]))
    ax = sns.stripplot(data=bond_lengths_inter, x="bond length (Å)", y="location", hue="design", order=order,
                       hue_order=design_order, orient="h", size=2, dodge=True)
    ax.set_aspect(0.5)
    ax.set_yticks(np.array(ax.get_yticks()) + 0.5, minor=True)
    ax.grid(axis="y", which="minor")
    plt.legend(bbox_to_anchor=(1.01, 1.0))
    plt.title("Planar length of O3'-P bond vs. location of base across origami designs")
    plt.savefig("bond_length_planar_vs_location.png", dpi=200)
    plt.close()


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
    COLORS = ["#ccbbff", 0.01, "yellow", 0.03, "red"]
    def color_map(value):
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
    SCAF_XO_COLOR = "#000000"
    STAP_XO_COLOR = "#c0c0c0"
    SCAF_XO = f"{terms.SCAF}_{terms.SCAF_XO}"
    STAP_XO = f"{terms.SCAF}_{terms.STAP_XO}"
    SCAF_TM = f"{terms.SCAF}_{terms.SCAF_TM}"
    STAP_TM = f"{terms.SCAF}_{terms.STAP_TM}"
    zorder = {"feature": 1, "base": 2, "text": 3}
    seqs = dict(zip(base_info["CanDo number"], base_info["base"]))
    locs = dict(zip(base_info["CanDo number"], base_info["location"]))
    xs = dict()
    ys = dict()
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    plt.axis("off")  # hide axes and ticks
    for edge, (double_helix1, double_helix2) in enumerate(edges):
        assert len(double_helix1) == len(double_helix2)
        for pos, (num1, num2) in enumerate(zip(double_helix1, double_helix2)):
            x = pos * X_INC
            y1 = -edge * Y_INC
            y2 = y1 - Y_SEP
            xs[num1], ys[num1] = x, y1
            xs[num2], ys[num2] = x, y2
            plt.scatter([x], [y1], c=[color_map(signals.get(num1, np.nan))], s=BASE_SIZE, zorder=zorder["base"])
            plt.scatter([x], [y2], c=[color_map(signals.get(num2, np.nan))], s=BASE_SIZE, zorder=zorder["base"])
            plt.text(x, y1, seqs[num1], zorder=zorder["text"], **SEQ_PARAMS)
            plt.text(x, y2, seqs[num2], zorder=zorder["text"], **SEQ_PARAMS)
            if num1 == 1 or num1 % NUM_PERIOD == 0:
                plt.text(x, y1 + NUM_SEP, num1, **NUM_PARAMS)
            if num2 == 1 or num2 % NUM_PERIOD == 0:
                plt.text(x, y2 - NUM_SEP, num2, **NUM_PARAMS)
            loc1 = locs[num1]
            loc2 = locs[num2]
            if loc1.startswith(SCAF_XO):
                # scaffold crossovers
                if loc1.endswith("5"):
                    # x position of the base to which it is crossing over
                    x_xo = xs.get(g_dn[num1])
                elif loc1.endswith("3"):
                    x_xo = xs.get(g_up[num1])
                else:
                    raise ValueError()
                if x_xo is not None:
                    plt.plot([x, x_xo], [y1, y2], c=SCAF_XO_COLOR, zorder=zorder["feature"])
            if loc2.startswith(SCAF_XO):
                # scaffold crossovers
                if loc2.endswith("5"):
                    # x position of the base to which it is crossing over
                    x_xo = xs.get(g_dn[num2])
                elif loc2.endswith("3"):
                    x_xo = xs.get(g_up[num2])
                else:
                    raise ValueError()
                if x_xo is not None:
                    plt.plot([x, x_xo], [y1, y2], c=SCAF_XO_COLOR, zorder=zorder["feature"])
            if loc1.startswith(STAP_XO):
                # staple crossovers
                if loc1.endswith("5"):
                    # x position of the base to which it is crossing over
                    x_xo = xs.get(g_ax[g_up[g_ax[num1]]])
                elif loc1.endswith("3"):
                    # x position of the base to which it is crossing over
                    x_xo = xs.get(g_ax[g_dn[g_ax[num1]]])
                else:
                    raise ValueError()
                if x_xo is not None:
                    plt.plot([x, x_xo], [y1, y2], c=STAP_XO_COLOR, zorder=zorder["feature"])
            if loc2.startswith(STAP_XO):
                # staple crossovers
                if loc2.endswith("5"):
                    x_xo = xs.get(g_ax[g_up[g_ax[num2]]])
                elif loc2.endswith("3"):
                    x_xo = xs.get(g_ax[g_dn[g_ax[num2]]])
                else:
                    raise ValueError()
                if x_xo is not None:
                    plt.plot([x, x_xo], [y1, y2], c=STAP_XO_COLOR, zorder=zorder["feature"])
            if loc1.startswith(SCAF_TM):
                # staple crossovers
                if loc1.endswith("5"):
                    x_tm = x - X_INC / 2
                elif loc1.endswith("3"):
                    x_tm = x + X_INC / 2
                else:
                    raise ValueError()
                if x_tm is not None:
                    plt.plot([x_tm, x_tm], [y1, y1 + TM_HEIGHT], c=SCAF_XO_COLOR, zorder=zorder["feature"])
            if loc2.startswith(SCAF_TM):
                # staple crossovers
                if loc2.endswith("5"):
                    x_tm = x + X_INC / 2
                elif loc2.endswith("3"):
                    x_tm = x - X_INC / 2
                else:
                    raise ValueError()
                if x_tm is not None:
                    plt.plot([x_tm, x_tm], [y2, y2 - TM_HEIGHT], c=SCAF_XO_COLOR, zorder=zorder["feature"])
            if loc1.startswith(STAP_TM):
                # staple crossovers
                if loc1.endswith("5"):
                    x_tm = x + X_INC / 2
                elif loc1.endswith("3"):
                    x_tm = x - X_INC / 2
                else:
                    raise ValueError()
                if x_tm is not None:
                    plt.plot([x_tm, x_tm], [y1, y1 + TM_HEIGHT], c=STAP_XO_COLOR, zorder=zorder["feature"])
            if loc2.startswith(STAP_TM):
                # staple crossovers
                if loc2.endswith("5"):
                    x_tm = x - X_INC / 2
                elif loc2.endswith("3"):
                    x_tm = x + X_INC / 2
                else:
                    raise ValueError()
                if x_tm is not None:
                    plt.plot([x_tm, x_tm], [y2, y2 - TM_HEIGHT], c=STAP_XO_COLOR, zorder=zorder["feature"])
    plt.savefig(fname, dpi=600)
