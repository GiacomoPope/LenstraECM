from gmpy2 import isqrt, ceil, mpz
from primesieve import primes

def trial_division(N):
    """
    Trial division, algorithm 8.1.1 (Cohen)
    """
    factors = {}
    # Trivial cases
    if N == 1:
        return {}, 1
    elif N == 2:
        return {2: 1}, 1
    elif N == 3:
        return {3: 1}, 1
    elif N == 4:
        return {2: 2}, 1
    elif N == 4:
        return {5: 1}, 1
    # initialise
    # Cohen suggests 500,000 but I think
    # Computers are faster now and we can
    # afford to go higher... Maybe even
    # Higher still...
    B = min(N, 1_000_000)
    P = primes(B)
    l = ceil(isqrt(N))
    # Trial divide
    for d in P:
        r = N % d
        while r == 0:
            N = N // d
            r = N % d
            l = ceil(isqrt(N))
            if d in factors:
                factors[d] += 1
            else:
                factors[d] = 1
            if N == 1:
                return factors, N
    if d >= l:
        factors[int(N)] = 1
        return factors, 1
    # Larger factors remain
    return factors, mpz(N)


if __name__ == '__main__':
    from gmpy2 import is_prime
    import random
    from factor_util import *

    for _ in range(1000):
        n = random.randint(0, 10**10)
        factors, hmm = trial_division(n)
        if not test_factors(n, factors):
            assert n % hmm == 0
            assert test_factors(n // hmm, factors)
