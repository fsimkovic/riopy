def count_res_in_frag(residues, min_frag_len=3):
    """Count the number of fragments in a list
    
    Parameters
    ----------
    residues : list
       A list of residue indexes
    min_frag_len : int
       The minimum length in residues of a fragment [default: 3]

    Returns
    -------
    int
       The number of fragments

    """

    if len(residues) < 1:
        return 0

    def valid_frag_len(c):
        return c >= min_frag_len

    res_count = 0
    frag_len = 1

    previous = residues[0]
    for current in residues[1:]:
        if current - 1 == previous:
            frag_len += 1
        else:
            if valid_frag_len(frag_len):
                res_count += frag_len
            frag_len = 1
        previous = current

    if valid_frag_len(frag_len):
        res_count += frag_len

    return res_count
