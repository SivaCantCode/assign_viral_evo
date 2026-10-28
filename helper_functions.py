def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    import math

    nucleotides = set("ACGT")
    length = 0
    mismatches = 0

    # skip indels (and ambiguous bases like N), compare the rest
    for ref_nt, dist_nt in zip(reference_sequence.upper(), distant_sequence.upper()):
        if ref_nt not in nucleotides or dist_nt not in nucleotides:
            continue
        length += 1
        if ref_nt != dist_nt:
            mismatches += 1

    if length == 0:
        return 0.0

    p = mismatches / length
    if p >= 0.75:
        return float("inf")

    # JC distance is per site, scale back to the number of sites
    d_per_site = -3 / 4 * math.log(1 - 4 / 3 * p)
    return d_per_site * length


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
