#!/usr/bin/env python3

import json
import argparse
import sys
class smith_waterman:
    """
    Smith Waterman alignment class

    ...

    Attributes
    -----------
    _seq1 : str
        A string containing the first sequence to be aligned.
    _seq2 : str
        A string containing the second sequence to be aligned.
    _match : double
        A positive real number that describes the gain in score
        when two bases match when performing alignment.
    _mismatch : double
        A negative real number that described the loss in score
        when two bases mismatch when performing alignment.
    _k_weight : double
        A positive real number that scales the length of a gap,
        deriving the penalty to the score when inserting a gap
        of arbitrary length during alignment.
    _solutions : list
        A list containing all the possible alignments.
    _mat : list
        The matrix containing the score of the Smith-Waterman
        alignment and a list of all the possible coordinates
        of the cell of origin.

    Methods
    -------
    align()
        Performs the Smith-Waterman alignment populating self._mat.
    solve(min_score = 0.0, n_sol = None, only_max = False, include_div = False print_format = "txt")
        Performs the traceback procedure for the solution with the top
        n_sol score with score >= min_score, or only the solutions
        with the top score if only_max is set to True and prints them
        in format of print_format. Divergent solutions are printed only
        when include_div = True

    """

    def __init__(self, seq1, seq2, match, mismatch, k_weight):
        """
        Parameters
        ----------
        seq1 : str
            The first sequence to align.
        seq2 : str
            The second sequence to align.
        match : double
            Gain in score for base match.
        mismatch : double
            Lost in score for base mismatch.
        k_weight : double
            Scale of lost in score for gap insertions.
        """

        self._seq_col = seq1
        self._seq_row = seq2
        self._match = match
        self._mismatch = mismatch
        self._k_weight = k_weight
        self._solutions = []
        self._mat = [[[0, []] for i in range(len(self._seq_row) + 1)] for j in range(len(self._seq_col) + 1)]

    def __score_matrix(self):
        """Prints the matrix containing the scores"""

        res = "\t\t" + "\t".join(self._seq_row) + "\n"
        for i in range(len(self._mat)):
            if i == 0:
                res += "\t" + "\t".join([str(j[0]) for j in self._mat[i]]) + "\n"
            else:
                res += self._seq_col[i - 1] + "\t" + "\t".join([str(j[0]) for j in self._mat[i]]) + "\n"
        return res

    def __traceback_matrix(self):
        """Prints the matrix containing the origins for cell's score"""

        res = "\t\t" + "\t".join(self._seq_row) + "\n"
        for i in range(len(self._mat)):
            if i == 0:
                res += "\t" + "\t".join([str(j[1:]) for j in self._mat[i]]) + "\n"
            else:
                res += self._seq_col[i - 1] + "\t" + "\t".join([str(j[1:]) for j in self._mat[i]]) + "\n"
        return res

    def __str__(self):
        """Print the smith-waterman object state"""

        res = "Match weight: {}\nMismatch weight: {}\nW_k: {}\n".format(self._match, self._mismatch, self._k_weight)
        res += "First sequence:\n{}\nSecond sequence:\n{}\n".format(self._seq_row, self._seq_col)
        res += "Matrix dimension: {}x{}".format(len(self._mat), len(self._mat[0])) + "\n"
        res += self.__score_matrix()
        res += self.__traceback_matrix()

        return(res)

    def __diagonal(self, up_value, base1, base2):
        """Computes the score of inserting a match or mismatch
        Parameters
        ----------
        up_value : double
            Score of the originating cell
        base1 : str
            Base in position in sequence 1
        base2 : str
            Base in position in sequence 2
        Return
        ------
        res : list
            List containing score and coordinates of origin
        """

        return up_value + (self._match if base1 == base2 else self._mismatch)

    def __gap_insert(self, col_values):
        """Compute the best scoring gap to insert when performing alignment
        Parameters
        ----------
        col_values : list
            Score list on where to insert a gap
        Return
        ------
        res : list
            List containing score and changed coordinate of gap origin
        """

        res = [-1, -1]
        for k in range(len(col_values)):
            if col_values[k] - self._k_weight*(len(col_values) - k) > res[0]:
                res = [col_values[k] - self._k_weight*(len(col_values) - k), k]

        return res

    def align(self):
        """Performs the alignment populating self._mat"""

        for i in range(1, len(self._mat)):
            for j in range(1, len(self._mat[i])):
                diagonal = [self.__diagonal(self._mat[i-1][j-1][0], self._seq_col[i-1], self._seq_row[j-1]), i - 1, j - 1]
                vert_gap = self.__gap_insert([self._mat[k][j][0] for k in range(i)]) + [j] # Adding coordinate that didn't change in gap insertion
                hor_gap = self.__gap_insert([self._mat[i][k][0] for k in range(j)])
                hor_gap = [hor_gap[0], i, hor_gap[1]] # Adding coordinate that didn't change in gap insertion
                # Populate matrix with best score between 0 (default) and one of the three operations
                max_score = max(diagonal[0], vert_gap[0], hor_gap[0], 0)
                if max_score > 0:
                    self._mat[i][j][0] = max_score
                    if diagonal[0] == max_score: # Base match or mismatch
                        self._mat[i][j][1].append([diagonal[1], diagonal[2]])
                    if vert_gap[0] == max_score: # Vertical gap
                        self._mat[i][j][1].append([vert_gap[1], vert_gap[2]])
                    if hor_gap[0] == max_score: # Horizontal gap
                        self._mat[i][j][1].append([hor_gap[1], hor_gap[2]])

    def __sort_end_point(self, end_point):
        """Function to sort end point by score"""

        return self._mat[end_point[0]][end_point[1]]

    def __get_all_solutions_end_points(self, min_score = 0.0, n_sol = None, only_max = False):
        """Returns all ending point for the top n_sol solution with score >= min_score
        Parameters
        ----------
        min_score = 0.0 : double
            Minimum score for a solution to be included.
        n_sol = None : int
            Maximum number of top solution to be returned.
        only_max = False : bool
            If true return only the solutions with max score.
        Return
        ------
        res : list
            List of all the possible ending point for solutions that respect filters.
        """
        if only_max:
            min_score = max(map(max, self._mat))[0]

        res = []
        for i in range(len(self._mat)):
            for j in range(len(self._mat[i])):
                if self._mat[i][j][0] >= min_score and self._mat[i][j][0] > 0:
                    res.append([i, j])

        res.sort(key = self.__sort_end_point, reverse = True) # Sort in descending score order
        if n_sol:
            res = res[:n_sol]
        return res

    def __single_sol_traceback(self, end_point, include_div = False, print_format = "txt"):
        """Performs the traceback procedure for a single endpoint
        Parameters
        ----------
        end_point : list
            List of all the possible end_points.
        Return
        ------
        res : list
            List containing all the resulting alignments.
            Each alignment is a dictionary containing different
            alignment qualities.
        """

        res = {"alignment_seq1" : "",
               "alignment_seq2" : "",
               "alignment_length" : 0,
               "score" : self._mat[end_point[0]][end_point[1]][0],
               "n_match" : 0,
               "n_mismatch" : 0,
               "n_gaps" : 0}

        end_point = [end_point]
        results = [res]
        t = 0
        while len(results) > 0: # t refers to results[t] originated from end_point[t]
            while self._mat[end_point[t][0]][end_point[t][1]][0] > 0: # Traceback for a single procedure
                local_res = results[t]
                temp_res = local_res.copy()
                div_number = 1
                if include_div: #Include divergent solution, if false is set to one as to take only the first origin for a cell score
                    div_number = len(self._mat[end_point[t][0]][end_point[t][1]][1])
                for i in range(div_number): # If the traceback diverges, add a solution to the result list
                    local_end_point = self._mat[end_point[t][0]][end_point[t][1]][1][i]
                    if i > 0:
                        local_res = temp_res.copy()
                        results.append(local_res)
                    # If score is moving up and on the left by one (base match or mismatch)
                    if local_end_point == [end_point[t][0] - 1, end_point[t][1] - 1]:
                        local_res["alignment_seq1"] = self._seq_row[local_end_point[1]] + local_res["alignment_seq1"]
                        local_res["alignment_seq2"] = self._seq_col[local_end_point[0]] + local_res["alignment_seq2"]
                        if self._seq_row[local_end_point[1]] == self._seq_col[local_end_point[0]]:
                            local_res["n_match"] += 1
                        else:
                            local_res["n_mismatch"] += 1
                    # If score is moving up (gap insertion)
                    elif local_end_point[0] == end_point[t][0] and local_end_point[1] != end_point[t][1]:
                        local_res["alignment_seq1"] = self._seq_row[end_point[t][1] - 1] + local_res["alignment_seq1"]
                        n_gap = end_point[t][1] - local_end_point[1]
                        local_res["alignment_seq2"] = "".join(["-" for k in range(n_gap)]) + local_res["alignment_seq2"]
                        local_res["n_gaps"] += n_gap
                    # If score is moving left (gap insertion)
                    elif local_end_point[0] != end_point[t][0] and local_end_point[1] == end_point[t][1]:
                        local_res["alignment_seq2"] = self._seq_col[end_point[t][0] - 1] + local_res["alignment_seq2"]
                        n_gap = end_point[t][0] - local_end_point[0]
                        local_res["alignment_seq1"] = "".join(["-" for k in range(n_gap)]) + local_res["alignment_seq1"]
                        local_res["n_gaps"] += n_gap
                    local_res["alignment_length"] += 1

                end_point[t] = self._mat[end_point[t][0]][end_point[t][1]][1][0] # Endpoint traceback
                if len(self._mat[end_point[t][0]][end_point[t][1]][1]) > 1:
                    end_point += self._mat[end_point[t][0]][end_point[t][1]][1][1:] # Add endpoint of divergent tracebacks

            #Print result to free memory for large alignments
            if print_format == "txt":
                print(self.__solution_to_string(results[t]))
            elif print_format == "tsv":
                print(self.__solution_to_tsv(results[t]))
            elif print_format == "json":
                print(self.__solution_to_json(results[t]))
            results.pop(0)
            end_point.pop(0)



        return results

    def solve(self, min_score = 0.0, n_sol = None, only_max = False, include_div = False, print_format = "txt"):
        """Performs the traceback procedure returning alignments
           with score >= min_score.
        Parameters
        ----------
        min_score = 0.0 : double
            The minimum score for a solution to be included.
        n_sol = None : int
            Compute the alignment only for the solution with the best n_sol scores.
        only_max = False : bool
            If true get only the solutions with the max score.
        include_div = False : bool
            If true print all divergent solutions
        print_format = "txt" : str
            Format in which solutions are printed
        """

        self.__print_header(print_format)

        end_points = self.__get_all_solutions_end_points(min_score, n_sol, only_max)
        for i in end_points:
            self._solutions += self.__single_sol_traceback(i, include_div, print_format)

    def __print_header(self, print_format):
        """Prints header according to print_format
        Parameters
        ----------
        print_format : str
            The format in which to print output header
        """

        if print_format == "txt":
            print("First sequence: {}\nSecond sequence: {}".format(self._seq_col, self._seq_row))
        elif print_format == "tsv":
            print("##First sequence: {}\n##Second sequence: {}".format(self._seq_col, self._seq_row))
            print("alignment_seq_1\talignemnt_seq2\talignment_length\tscore\tmatch\tmismatch\tgaps")
        elif print_format == "json":
            print(json.dumps({"First_sequence" : self._seq_col}, indent = 2))
            print(json.dumps({"Second_sequence" : self._seq_row}, indent = 2))

    def __solution_to_json(self, sol):
        """Converts solution into a json object.
        Return
        ------
        res : json
            Solution converted into a json object.
        """

        res = json.dumps(sol, indent = 2)

        return res

    def __solution_to_tsv(self, sol):
        """Converts solution into a tsv-like row.
        Return
        ------
        res : json
            Solution converted into a tsv-like string.
        """
        res = '\t'.join(str(sol[val]) for val in sol)
        return res

    def __solution_to_string(self, sol):
        """Transforms a solution in a printable string.
        Parameters
        ----------
        sol : dict
            The solution to be transformed into string.
        Return
        ------
        res : str
            The solution as a string.
        """

        res = "Alignment:\n"
        res += "\t" + sol["alignment_seq1"]
        res += "\n\t" + sol["alignment_seq2"] + "\n"
        res += "Score: " + str(sol["score"]) + "\n"
        res += "Alignment length: " + str(sol["alignment_length"]) + "\n"
        res += "Number of matches: " + str(sol["n_match"]) + "\n"
        res += "Number of mismatches: " + str(sol["n_mismatch"]) + "\n"
        res += "Number of gaps: " + str(sol["n_gaps"]) + "\n"

        return res

if __name__ == "__main__":
    # Input management
    parser = argparse.ArgumentParser()
    parser.add_argument("first_sequence", help = "The reference sequence")
    parser.add_argument("second_sequence", help = "The sequence to be aligned")
    parser.add_argument("-m", "--match", help = "Score of a base match [default 3]", type = float, default = 3)
    parser.add_argument("-p", "--penalty", help = "Penalty of a base match [default -3]", type = float, default = -3)
    parser.add_argument("-g", "--gap", help = "Scale for the penalty of a gap insertion [default 2]", type = float, default = 2)
    parser.add_argument("-s", "--min-score", help = "Minimum score for a solution to be included [default 0]", type = float, default = 0)
    parser.add_argument("-M", "--only-max", help = "Print only the max scoring solutions [default False]", action = "store_true")
    parser.add_argument("-d", "--include-divergent", help = "Print all the divergent solution [default False]", action = "store_true")
    parser.add_argument("-n", "--n-result", help = "Print only the alignment with the best N scores [default all]", type = int)
    parser.add_argument("-f", "--format", help = "Specify output format [default txt]", type = str, default = "txt", choices = ["txt", "json", "tsv"])
    parser.add_argument("-o", "--output-file", help = "Specify output file [default stdout]", type = str)
    args = parser.parse_args()
    # End of input management

    # Redirect output if needed
    if args.output_file:
        sys.stdout = open(args.output_file, 'w')

    # Initialize alignment object
    alignment_out = smith_waterman(args.first_sequence, args.second_sequence, args.match, args.penalty, args.gap)
    # Populate score matrix
    alignment_out.align()
    # Perform traceback procedure and print results
    alignment_out.solve(args.min_score, args.n_result, args.only_max, args.include_divergent, args.format)
