def test_factors(N, factors):
    """
    Take all the factors and exponents and make
    sure nothing funky has gone on during factoring...
    """
    test = 1
    for p, e in factors.items():
        test *= p**e
    return N == test
