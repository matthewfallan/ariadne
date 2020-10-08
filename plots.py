"""
ARIADNE - Plots module

Make figures.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


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

