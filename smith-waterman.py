#!/usr/bin/env python3

import json
import argparse

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
    solve(min_score = 0.0, n_sol = None)
        Performs the traceback procedure for the solution with the top
        n_sol score with score >= min_score.
    solutions_to_json()
        Returns all the computed solutions in a json format.
    solutions_to_tsv()
        Returns all the computed solutions in a tsv format.
    get_solutions()
        Returns all the solutions
    solutions_to_string()
        Returns all the solutions in a printable string.
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

    def __get_all_solutions_end_points(self, min_score = 0.0, n_sol = None):
        """Returns all ending point for the top n_sol solution with score >= min_score
        Parameters
        ----------
        min_score = 0.0 : double
            Minimum score for a solution to be included.
        n_sol = None : int
            Maximum number of top solution to be returned.
        Return
        ------
        res : list
            List of all the possible ending point for solutions that respect filters.
        """

        res = []
        for i in range(len(self._mat)):
            for j in range(len(self._mat[i])):
                if self._mat[i][j][0] >= min_score and self._mat[i][j][0] > 0:
                    res.append([i, j])

        res.sort(key = self.__sort_end_point, reverse = True) # Sort in descending score order
        if n_sol:
            res = res[:n_sol]
        return res

    def __single_sol_traceback(self, end_point):
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

        res = {"seq1" : self._seq_row,
               "seq2" : self._seq_col,
               "alignment_seq1" : "",
               "alignment_seq2" : "",
               "alignment_length" : 0,
               "score" : self._mat[end_point[0]][end_point[1]][0],
               "n_match" : 0,
               "n_mismatch" : 0,
               "n_gaps" : 0}

        end_point = [end_point]
        results = [res]
        t = 0
        while t < len(end_point): # t refers to results[t] originated from end_point[t]
            while self._mat[end_point[t][0]][end_point[t][1]][0] > 0: # Traceback for a single procedure
                local_res = results[t]
                temp_res = local_res.copy()
                for i in range(len(self._mat[end_point[t][0]][end_point[t][1]][1])): # If the traceback diverges, add a solution to the result list
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

            t += 1

        return results

    def solve(self, min_score = 0.0, n_sol = None):
        """Performs the traceback procedure returning alignments
           with score >= min_score.
        Parameters
        ----------
        min_score = 0.0 : double
            The minimum score for a solution to be included.
        n_sol = None : int
            Compute the alignment only for the solution with the best n_sol scores.
        """

        end_points = self.__get_all_solutions_end_points(min_score, n_sol)
        for i in end_points:
            self._solutions += self.__single_sol_traceback(i)

    def solutions_to_json(self):
        """Converts list of solution into a json object.
        Return
        ------
        res : json
            Solution list converted into a json object.
        """

        res = json.dumps(self._solutions, indent = 2)

        return res

    def solutions_to_tsv(self):
        """Converts list of solution into a tsv-like string.
        Return
        ------
        res : json
            Solution list converted into a tsv-like string.
        """
        keys = self._solutions[0].keys()

        result = [list(keys)] + [list(row.values()) for row in self._solutions]
        res = '\n'.join(['\t'.join(str(val) for val in row) for row in result])
        return res

    def get_solutions(self):
        """Returns all the solutions
        Return
        ------
        self._solutions : list
            returns all the solutions
        """

        return self._solutions

    def __solution_to_string(self, sol):
        """Transforms a solution in a printable string
        Parameters
        ----------
        sol : dict
            The solution to be transformed into string
        Return
        ------
        res : str
            The solution as a string
        """

        res = "First sequence: " + sol["seq1"] + "\n"
        res += "Second sequence: " + sol["seq2"] + "\n"
        res += "Alignment:\n"
        res += "\t" + sol["alignment_seq1"]
        res += "\n\t" + sol["alignment_seq2"] + "\n"
        res += "Score: " + str(sol["score"]) + "\n"
        res += "Alignment length: " + str(sol["alignment_length"]) + "\n"
        res += "Number of matches: " + str(sol["n_match"]) + "\n"
        res += "Number of mismatches: " + str(sol["n_mismatch"]) + "\n"
        res += "Number of gaps: " + str(sol["n_gaps"]) + "\n"

        return res

    def solutions_to_string(self):
        """Return all solutions in text format
        Return
        ------
        res : str
            All the solutions in text format
        """
        res = ""
        for i in self._solutions:
            res += self.__solution_to_string(i)

        return res


if __name__ == "__main__":
    # Input management
    parser = argparse.ArgumentParser()
    parser.add_argument("first_sequence", help = "The first sequence to be aligned")
    parser.add_argument("second_sequence", help = "The second sequence to be aligned")
    parser.add_argument("-m", "--match", help = "Score of a base match [default 3]", type = float, default = 3)
    parser.add_argument("-p", "--penalty", help = "Penalty of a base match [default -3]", type = float, default = -3)
    parser.add_argument("-g", "--gap", help = "Scale for the penalty of a gap insertion [default 2]", type = float, default = 2)
    parser.add_argument("-s", "--min-score", help = "Minimum score for a solution to be included [default 0]", type = float, default = 0)
    parser.add_argument("-n", "--n-result", help = "Print only the alignment with the best N score [default all]", type = int)
    parser.add_argument("-f", "--format", help = "Specify output format [default txt]", type = str, default = "txt", choices = ["txt", "json", "tsv"])
    parser.add_argument("-o", "--output-file", help = "Specify output file [default stdout]", type = str)
    args = parser.parse_args()
    # End of input management

    # Initialize alignment object
    alignment_out = smith_waterman(args.first_sequence, args.second_sequence, args.match, args.penalty, args.gap)
    # Populate score matrix
    alignment_out.align()
    # Perform traceback procedure
    alignment_out.solve(args.min_score, args.n_result)


    # Output format selections
    out = ""
    if args.format == "txt":
        out = alignment_out.solutions_to_string()
    elif args.format == "json":
        out = alignment_out.solutions_to_json()
    elif args.format == "tsv":
        out = alignment_out.solutions_to_tsv()
    print(alignment_out.get_solutions())
    # Output redirection (if needed)
    if not args.output_file:
        print(out)
    else:
        with open(args.output_file, "w") as out_file:
            out_file.write(out)
