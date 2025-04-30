#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=iecor-qtpie
#SBATCH --ntasks-per-node=1
#SBATCH --partition=eng-instruction
#SBATCH --mem=512000
#SBATCH --account=25sp-cs581a-eng

# Default value for quartet
quartet=11
guidetree=""
name="tmp"
input=""

# Function to show usage
usage() {
  echo "Usage: $0 -g <guidetree> -i <input> [-n <name>] [-q <quartet>]"
  echo "  -g, --guidetree   Path to the guide tree file (required)"
  echo "  -i, --input       Path to the input file (required)"
  echo "  -q, --quartet     Quartet value (optional, default: 11)"
  echo "  -n, --name        Name of the folder to put results in (optional, default: tmp)"
  exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -g|--guidetree)
      guidetree="$2"
      shift 2
      ;;
    -i|--input)
      input="$2"
      shift 2
      ;;
    -q|--quartet)
      quartet="$2"
      shift 2
      ;;
    -n|--name)
      name="$2"
      shift 2
      ;;
    -*)
      echo "Unknown option: $1"
      usage
      ;;
    *)
      echo "Unexpected argument: $1"
      usage
      ;;
  esac
done

# Check for required arguments
if [[ -z "$guidetree" || -z "$input" ]]; then
  echo "Error: --guidetree and --input are required"
  usage
fi


HASH=$(uuidgen)

QUARTET_PATH=~/scratch/qtpie_tmpqs_$HASH.nwk
time python3 $TALLIS/QuartetMethods/scripts/py/printQuartets.py\
    -i $input\
    -q 12\
    > $QUARTET_PATH

echo Quartet Generation on EVANS-ONE complete. $(cat $QUARTET_PATH | wc -l) quartets.

time python3 $TALLIS/QTPIE/runQTPIE.py\
    -g $guidetree\
    -q $QUARTET_PATH\
    -A /projects/illinois/eng/cs/warnow/zxliu2/OneMostProb/ASTRAL/astral.5.7.8.jar\
    -o $TALLIS/QTPIE/$name\
    > $TALLIS/QTPIE/$name/qtpie_tree.nwk 2> $TALLIS/QTPIE/$name/qtpie_tree.log
     