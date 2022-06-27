# Smith Waterman implementation
A Smith Waterman implementation in python, done for the course of algorithms for bioinformatics, of professor Blanzieri, QCB 2021-2022.


## How to use
After cloning the repo to perform an alignment between sequences:

```smith-waterman.py [-h] [-m MATCH] [-p PENALTY] [-g GAP] [-s MIN_SCORE] [-M] [-n N_RESULT] [-f {txt,json,tsv}] [-o OUTPUT_FILE] first_sequence second_sequence```

### Options:

- First sequence [required] the reference to which align the second sequence.
- Second sequence [required] the sequence to align to the first one.
- -h, --help [optional] print the help message and exit.
- -m, --match [optional] change the score for a base match [default 3]
- -p, --penalty [optional] change the penalty for a base match [default -3]
- -g, --gap [optional] scale for the penalty for a gap insertion [default 2]
- -s, --min-score [optional] the minumum score for a solution to be included [default 0]
- -M, --only-max [optional] print only the max scoring solutions [default False]
- -d, --include-divergent [optional] print divergent solutions [default True] (WARNING this solution is exponential depending on the length of the alignment, use it at your own risk.
- -n, --n-results [optional] print only the alignment with the best N scores [default all]
- -f, --format [optional] choose the output format between tsv, json or txt [default txt]
- -o, --output-file [optional] specify the file in which to save the alignemtn [default stdout]
