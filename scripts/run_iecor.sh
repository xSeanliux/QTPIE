#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=iecor-qtpie
#SBATCH --ntasks-per-node=1
#SBATCH --partition=eng-instruction
#SBATCH --mem=512000
#SBATCH --account=25sp-cs581a-eng

QUARTET_PATH=~/scratch/quartets_iecor.nwk
time python3 $TALLIS/QuartetMethods/scripts/py/printQuartets.py\
    -i /projects/illinois/eng/cs/warnow/zxliu2/QuartetMethods/example/iecor/iecor_screened_question_marks.csv\
    -q 12\
    > $QUARTET_PATH

echo Quartet Generation on EVANS-ONE complete. $(cat $QUARTET_PATH | wc -l) quartets.

time python3 $TALLIS/QTPIE/runQTPIE.py\
    -g /projects/illinois/eng/cs/warnow/zxliu2/QTPIE/test/guide_iecor.nwk\
    -q $QUARTET_PATH\
    -A /projects/illinois/eng/cs/warnow/zxliu2/OneMostProb/ASTRAL/astral.5.7.8.jar\
    -o $TALLIS/QTPIE/tmp\
    > $TALLIS/QTPIE/tmp/IECOR-quartets.nwk 2> $TALLIS/QTPIE/tmp/IECOR-quartets.log
     