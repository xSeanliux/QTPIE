# Quartet Trees Polytomy InferencE (QTPIE)
Resolving polytomies with quartets. Can be used to scale up ASTRAL (or other quartet methods) with a guide tree.

## Usage 
```bash
python3 runQuartetParsimony.py\
    -g [GUIDETREE].nwk\
    -f [FORMAT]\
    -q [QUARTETS].nwk\
    -A [ASTRALPATH]\
    -o [OUTPUT]\
    2> run_log.log > output.nwk
```
Where `[GUIDETREE]` is the path to your guide tree, which is in format `[FORMAT]` (defaults to `newick`, for valid values see the [BioPython documentation](https://biopython.org/wiki/Phylo)); `[QUARTETS]` is a list of quartets in newick format. Each taxon name must be alphanumeric (matching `\w+`), and be of the format `((a,b),(c,d));` (see `test/` for quartet list examples). The quartet parsing algorithm uses regular expressions so be careful! The `[ASTRALPATH]` is the path the the [ASTRAL-III](https://github.com/smirarab/ASTRAL) `.jar` file. Finally, `[OUTPUT]` is the folder where (temporary) outputs will be. Best to make it a temporary folder.