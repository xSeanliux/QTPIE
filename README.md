# Quartet Trees Polytomy InferencE (QTPIE)
Resolving polytomies with quartets. Can be used to scale up ASTRAL (or other quartet methods) with a guide tree.

## Setup 
### Python dependencies 
This project requires Python. It was developed on Python 3.12.7, and depends on the following packages: 
- [treeswift](https://github.com/niemasd/TreeSwift) (v1.1.45)
- [BioPython](https://biopython.org/) (v1.78)
### ASTRAL 
Getting ASTRAL is easy, as it is already a submodule within QTPIE! Simply run 
```bash
git submodule init; git submodule update
```
to download the files in the `ASTRAL/` folder, and then build it by doing 
```bash
cd ASTRAL/; bash ./make.sh; cd ..
```
Verify that ASTRAL is installed with the following command: 
```bash
[[ -f ASTRAL/astral.5.7.8.jar ]] && echo "ASTRAL properly installed!"
```

## Usage 
### (UNDER CONSTRUCTION, DOES NOT WORK YET) I have a CSV file of linguistic data
Do you have a set of data like in `data/rt_2025_poly_screened_lv_1.csv` and a guide tree like in `data/guide_rt.nwk`? Great! Use the bash script in `run_qtpie.sh`. 
- Author's note: does not work yet because quartet generation scheme is not finalised. Please stay tuned! In the meantime feel free to go over to [OneMostProb](https://github.com/xSeanliux/OneMostProb) to see how quartet generation will be done (WIP).
### I already have quartets
**NOTE**: Unlike ASTRAL, QTPIE does *not* accept gene trees. Only quartet trees are permitted as of the moment! In addition, quartet trees have to be inputted in a specific format (as I use regular expressions instead of a full newick parser to speed up quartet ingestion).
```bash
python3 runQuartetParsimony.py\
    -g [GUIDETREE].nwk\
    -f [FORMAT]\
    -q [QUARTETS].nwk\
    -A [ASTRALPATH]\
    -o [OUTPUT]\
    2> run_log.log > output.nwk
```
Where 
- `[GUIDETREE]` is the path to your guide tree, which is in format `[FORMAT]` (defaults to `newick`, for valid values see the [BioPython documentation](https://biopython.org/wiki/Phylo)); 
- `[QUARTETS]` is a list of quartets in newick format. Each taxon name must be alphanumeric (matching `\w+`), and be of the format `((a,b),(c,d));` (see `test/` for quartet list examples). The quartet parsing algorithm uses regular expressions so be careful! 
- Finally, `[OUTPUT]` is the folder where (temporary) outputs will be. Best to make it a temporary folder.
- (Optional!) The `[ASTRALPATH]` is the path the the [ASTRAL-III](https://github.com/smirarab/ASTRAL) `.jar` file. 