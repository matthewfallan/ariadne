# ariadne
Analyze outputs from DAEDALUS.

Usage:
python ariadne.py path/to/directory/of/DAEDALUS/outputs/for/one/design [path/to/directory/of/DAEDALUS/outputs/for/a/second/design] ...

Outputs:
Ariadne creates a subdirectory called "ariadne" within the directory containing the outputs from DAEDALUS.
Two files are output into that directory: base_info.tsv and chimerax_color.txt

base_info.txt
A tabular file listing every base in the scaffold and staples, the ID of each base in the CanDo and PDB files, the bond lengths in that base, and the location of each base in the origami. Locations are comma-delimited strings of the strand on which the base is located ("scaffold" or "staple"), the structural feature of that base ("staple/scaffold crossover", "staple/scaffold/edge terminus", "vertex", or "middle" that is not any of the other features), and which side of the feature the base is on ("5" for 5', "3" for 3', "0" if the side does not exist).

chimerax_color.txt
A script that can be run by Chimera X to open the PDB file and color-code each base based on the structural feature at that base. Staple crossovers are green, scaffold crossovers are teal, staple termini are pink, scaffold termini are brown, edge termini are blue, vertex bases are yellow, and middle bases (i.e. none of the above) are gray.
