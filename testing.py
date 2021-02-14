import sys
import numbers
import io
import recombinant
import submatrix

# The commonly-used BLOSUM62 substitution matrix for amino acid sequences
blosum62 = submatrix.read_substitution_matrix("BLOSUM62.txt")

# A simple match=+1 and mismatch=-1 substitution matrix for DNA
basic_dna_submatrix = submatrix.match_mismatch_matrix(1, -1, "ACGT")

# BLOSUM62 substitution matrix with space score = -11
blosum62_s11 = submatrix.submatrix_with_spaces(blosum62, -11)

# BLOSUM62 substitution matrix with space score = -1
blosum62_s1 = submatrix.submatrix_with_spaces(blosum62, -1)

# Simple DNA matrix with space score = -1
basic_dna_s1 = submatrix.submatrix_with_spaces(basic_dna_submatrix, -1)

# call the methods from recombine.py
def score_recombine_align(alignment, parents, substitution_matrix, r):
    return recombinant.score_recombine_align(alignment, parents, substitution_matrix, r)


def recombine_align(x, y, substitution_matrix, r):
    return recombinant.recombine_align(x, y, substitution_matrix, r)


# Dictionary mapping test case names to inputs
test_case_inputs = {
    'small_1':  ("MM",
                 ["MM",
                  "DD"],
                 blosum62_s11, -15),
    'small_2':  ("DD",
                 ["MM",
                  "DD"],
                 blosum62_s11, -15),
    'small_3':  ("DM",
                 ["MM",
                  "DD"],
                 blosum62_s11, -7),
    'small_4':  ("DM",
                 ["MM",
                  "DD"],
                 blosum62_s11, -9),
    'small_5':  ("DIDMM",
                 ["MMMM",
                  "DDDD"],
                 blosum62_s11, -7),
    'small_6':  ("DDMM",
                 ["MMM-M",
                  "DDDDD"],
                 blosum62_s11, -7),
    'small_7':  ("DDMM",
                 ["MMMM-M",
                  "D-DDDD"],
                 blosum62_s11, -7),
    'small_8':  ("MMDDMM",
                 ["MMMMMM",
                  "DDDDDD"],
                 blosum62_s11, -7),
    'small_9':  ("DMDMDM",
                 ["MMMMMM",
                  "DDDDDD"],
                 blosum62_s11, -8),
    'small_10': ("M",
                 ["MMMM",
                  "DDDD"],
                 blosum62_s11, -7),
    'small_11': ("WBW",
                 ["AAA",
                  "BBB"],
                 blosum62_s1, -7),
    'small_12': ("ABWB",
                 ["AAA",
                  "BBB"],
                 blosum62_s11, -8),
    'large_1':  ("QTRAGCLIGCLHVNSPTNSPRRARSVA",
                 ["QTRAGC-LIGAEHVNNSTNSRSVA",
                  "PS--GCSLSGCLH--SPTNSRSVS"],
                 blosum62_s11, -7),
    'large_2':  ("AAGGGTCGCACGGCTCCG",
                 ["AAG--GCAGGCGATCGG--TCCG",
                  "ACGTAGGTCGC-AACGGCATGCG"],
                 basic_dna_s1, -3)
}

# Dictionary mapping test case names to correct outputs
test_case_correct_outputs = {
    'small_1': (10,   ['MM',
                       'MM',
                       'DD'],
                '00'),
    'small_2': (12,   ['DD',
                       'MM',
                       'DD'],
                '11'),
    'small_3': (4,    ['DM',
                       'MM',
                       'DD'],
                '10'),
    'small_4': (3,    ['DM',
                       'MM',
                       'DD'],
                '11'),
    'small_5': (4,    ['DIDMM',
                       'M-MMM',
                       'D-DDD'],
                '11100'),
    'small_6': (15,   ['DDM-M',
                       'MMM-M',
                       'DDDDD'],
                '11000'),
    'small_7': (15,   ['D-DM-M',
                       'MMMM-M',
                       'D-DDDD'],
                '111000'),
    'small_8': (18,   ['MMDDMM',
                       'MMMMMM',
                       'DDDDDD'],
                '001100'),
    'small_9': (9,    ['DMDMDM',
                       'MMMMMM',
                       'DDDDDD'],
                '111110'),
    'small_10': (-28, ['---M',
                       'MMMM',
                       'DDDD'],
                 '0000'),
    'small_11': (0,   ['--WBW',
                       'AA-A-',
                       'BB-B-'],
                 '11111'),
    'small_12': (-5,   ['ABWB',
                        'AA-A',
                        'BB-B'],
                 '1111'),
    'large_1': (54, ["QTRAGC-LIGCLHVNSPTNSPRRARSVA",
                     "QTRAGC-LIGAEHVNNSTNS----RSVA",
                     "PS--GCSLSGCLH--SPTNS----RSVS"],
                "0000000001110001111111111111"),
    'large_2': (11, ["AAG--GGTCGC--ACGGC-TCCG",
                     "AAG--GCAGGCGATCGG--TCCG",
                     "ACGTAGGTCGC-AACGGCATGCG"],
                "00000111111111111100000")
}


def remove_gaps(s):
    """Returns a string with space characters ('-') removed"""
    return s.replace('-', '')


def transpose_alignment(alignment):
    """Returns a column-based alignment from a row-based alignment or vice versa"""
    return list(map(''.join, zip(*alignment)))


def remove_all_space_cols(alignment):
    """Returns an alignment with columns containing all space characters removed"""
    all_space_col = '-' * len(alignment)
    return transpose_alignment(col for col in transpose_alignment(alignment) if col != all_space_col)


def check_valid_alignment_result(result, x, y):
    """Checks that the alignment result is valid for sequences x and y."""
    assert isinstance(result, tuple), "Output is not a tuple"
    assert len(result) == 3, "Output does not have exactly three elements"

    score, alignment, parents = result
    assert isinstance(score, numbers.Number), "Score is not a number"
    assert isinstance(alignment, list), "Alignment is not a list"
    assert isinstance(parents, str), "Parents is not a string"

    assert len(alignment) == 3, "Alignment does not have exactly three elements"
    assert all(isinstance(element, str)
               for element in alignment), "Alignment elements are not strings"
    num_columns = len(alignment[0])
    assert all(len(
        row) == num_columns for row in alignment), "Alignment strings do not have the same length"

    assert remove_gaps(alignment[0]) == x, "First string of alignment is not x"
    assert remove_gaps(alignment[1]) == remove_gaps(
        y[0]), "Second string of alignment is not y_0"
    assert remove_gaps(alignment[2]) == remove_gaps(
        y[1]), "Third string of alignment is not y_1"
    assert remove_all_space_cols(
        alignment[1:]) == y, "Pairwise alignment of y_0 and y_1 is not preserved"

    assert len(
        parents) == num_columns, "Length of parent string is not equal to the number of alignment columns"
    assert set(parents).issubset(
        {'0', '1'}), "Parent string does not consist of only '0' and '1'"


def check_valid_alignment_score(result, submatrix, r):
    """Checks that the computed score of the alignment is equal to the score given in the result."""
    score, alignment, parents = result
    computed_score = score_recombine_align(alignment, parents, submatrix, r)
    assert computed_score == score, "Computed score ({}) does not equal the returned score ({})".format(
        computed_score, score)


def pprint_recombine_align(recombine_align_result, names=["x", "y0", "y1"], width=80, stream=sys.stdout):
    """Pretty prints the result of recombine_align.

    Args:
        recombine_align_result: an output from recombine_align, which is a tuple (score, alignment, parents)
        names: a list of three strings giving the names of the three sequences x, y_0, and y_1
        width: an integer specifying the number of columns of the alignment to print before wrapping
        stream: the stream (a file-like object) to which to print
    """
    score, alignment, parents = recombine_align_result
    num_columns = len(parents)
    name_width = max(map(len, names))
    all_names = names + [""]
    all_rows = alignment + [parents]

    print(f"Score: {score}", file=stream)
    for start_column in range(0, num_columns, width):
        for name, row in zip(all_names, all_rows):
            row_slice = row[start_column: start_column + width]
            print(f"{name:>{name_width}} {row_slice}", file=stream)
        print(file=stream)


def check_test_case(case_name, test_name=None, valid_result=True, valid_score=True, correct_alignment=True):
    inputs = test_case_inputs[case_name]
    correct_output = test_case_correct_outputs[case_name]
    result = recombine_align(*inputs)
    if valid_result:
        check_valid_alignment_result(result, *inputs[:2])
    if valid_score:
        check_valid_alignment_score(result, *inputs[2:])
    if correct_alignment:
        error_message = io.StringIO()
        print("Incorrect output", file=error_message)
        print("your function returned:", file=error_message)
        pprint_recombine_align(result, stream=error_message)
        print("the correct output is:", file=error_message)
        pprint_recombine_align(correct_output, stream=error_message)
        assert result == correct_output, error_message.getvalue()
    print("SUCCESS:", test_name if test_name else case_name, "passed!")


# TEST: small_1 valid output
check_test_case("small_1", test_name="small_1 valid output",
                valid_score=False, correct_alignment=False)

# TEST: small_1 valid score
check_test_case("small_1", test_name="small_1 valid score",
                correct_alignment=False)

# TEST: small_1 correct
check_test_case("small_1", test_name="small_1 valid score")

# TEST: small_2
check_test_case("small_2")

# TEST: small_3
check_test_case("small_3")

# TEST: small_4
check_test_case("small_4")

# TEST: small_5
check_test_case("small_5")

# TEST: small_6
check_test_case("small_6")

# TEST: small_7
check_test_case("small_7")

# TEST: small_8
check_test_case("small_8")

# TEST: small_9
check_test_case("small_9")

# TEST: small_10
check_test_case("small_10")

# TEST: small_11
check_test_case("small_11")

# TEST: small_12
check_test_case("small_12")

# TEST: large_1
check_test_case("large_1")

# TEST: large_2
check_test_case("large_2")

# TEST: large_3
# BEGIN HIDDEN TESTS
test_case_inputs['large_3'] = (
    'CGACGCATAGTGACACGCAAACTTGACCGG',
    ['TATCTCAGCCCCT-AAGATTGCGCAATAGC',
     '-C-T-GCTTGCGCTTC--TTCACACCTGGC'],
    basic_dna_s1, -2)

test_case_correct_outputs['large_3'] = (
    1,
    ['-CGACGCATAGTGACACGC-AA-ACTTG--ACC-GG-',
     'TATCTCA-G--CC---CCT-AAGA-TTGCGCAATAGC',
     '-C-T-GC-T--TG---CGCTTC---TTCACACCTGGC'],
    '1111111111111111111000000000111111111')
check_test_case("large_3")
# END HIDDEN TESTS

# TEST: large_4
# BEGIN HIDDEN TESTS
test_case_inputs['large_4'] = (
    'TGCCTATATGCGAATGAGGCGGTCTGGTAC',
    ['AACA-CTTTCG-TGGACGTGGGATATCG',
     'CT-CTCGGG-AACTT-ACATGA-A-T--'],
    basic_dna_s1,
    -2)

test_case_correct_outputs['large_4'] = (
    2,
    ['TGC--CTATATGCGAATGAGGCGGTCTGG-TA-C-',
     'AACA-CT-T-T-CG--TG-GAC-GT-GGGATATCG',
     'CT-CTCG-G-G--A-ACT-T-A-CA-TGA-A-T--'],
    '00000000000000000000000000000000000')

check_test_case("large_4")
# END HIDDEN TESTS

# TEST: large_5
# BEGIN HIDDEN TESTS
test_case_inputs['large_5'] = (
    'CCTCAACGGGGATCCACGTTCGATTACTTG',
    ['C-CCGCCAATC-AGACT-ATGCCGATTG-A',
     'CGACC-TG--GATGTA-CTTTCAGGGACTC'],
    basic_dna_s1, -2)

test_case_correct_outputs['large_5'] = (
    1,
    ['C-CTCAACGGGGATCCA-CGTTCGATTACTTG--',
     'C-C-CGCCAATC-AGACT-ATGC-C-GA-TTG-A',
     'CGA-CC-TG--GATGTA-CTTTC-A-GG-GACTC'],
    '0000000011111111111111111000000000')

check_test_case("large_5")
# END HIDDEN TESTS

# TEST: large_6
# BEGIN HIDDEN TESTS
test_case_inputs['large_6'] = (
    'ARYDHYAWPGKDTECWRERLHQLVHFLLVF',
    ['TMIGKDKFRVHGIMDHWQLHAGCWVCIRYQ',
     'MPK-YNYLWCIPRYTRRSREAFGWKMGYCE'],
    blosum62_s11,
    -2)

test_case_correct_outputs['large_6'] = (
    -13,
    ['ARY-DHYAWPGKDTECWRERLHQLVHFLLVF',
     'TMIGKDKFRVHGIMDHWQLHAGCWV-CIRYQ',
     'MPK-YNYLWCIPRYTRRSREAFGWK-MGYCE'],
    '0001111110000000001111100000000')

check_test_case("large_6")
# END HIDDEN TESTS

# TEST: large_7
# BEGIN HIDDEN TESTS
test_case_inputs['large_7'] = (
    'AMYIQMGQMWHYPNYWKAEMGAKLSFYCKVNAVYFTKGEYLGLAMQLCDWRWHDMIQDVQKKRISAVAAFNKCGHSHHMFIEKRPYEMCVTQSSCEQNYW',
    ['VSVHNFV-WTVNKMFYAACMKDMSHLYKYSH-ILMEFPTNWLPIQYDPKCPMLNWYLIKDWDNKRDSPKANEHGLMPIMRVFFT-SAFTILSVIQWAY-',
     'WLCNDCGICHMTASYSYTQLAHKYRTEMLSVYAGRQMLYNMFLLSASHAMDW-CDGMLVHDLNWKHMATKIVWQNHETWASEWCCSMHLRDFMHLTDDC'],
    blosum62_s11, -2)

test_case_correct_outputs['large_7'] = (
    22,
    ['AMYIQMGQMWHY-PNYWKAEMGAKLSFY-CKV-NAVYFTKGEYLGLAMQLCDW-RWHDMIQDVQKKRISAVAAFNKCGHSHHMFIEKRPYEMCVTQSSCEQNYW',
     'VSVHNFV--WTVNKMFYAACMKDMSHLYKYSH-ILMEFPTNWLPIQYDPKCPMLNWY-LIKDWDNKRDSPK-A-NEHG-LMPIMRVFFT-SAFTILSVIQWAY-',
     'WLCNDCG-ICHMTASYSYTQLAHKYRTEMLSVYAGRQMLYNMFLLSASHAMDW-CDG-MLVHDLNWKHMAT-K-IVWQ-NHETWASEWCCSMHLRDFMHLTDDC'],
    '01111111101111110001111110001111000000001111111111011100000000000000011000000011111111111111111100000001')

check_test_case("large_7")
# END HIDDEN TESTS
