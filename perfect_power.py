from primesieve import primes
from gmpy2 import iroot_rem

def perfect_power(n):
    """
    Assumes n is a perfect power
    I check this with gmpy's is_power
    """
    e_max = n.bit_length()
    for e in reversed(primes(e_max)):
        r, rem = iroot_rem(n, e)
        if not rem:
            _r, _e = perfect_power(r)
            if _r:
                r, e = _r, e*_e
            assert pow(r,e) == n
            return int(r), e
    return None, None


if __name__ == '__main__':
    from gmpy2 import mpz
    from time import time
    n = mpz(3**1024)
    t = time()
    print(perfect_power(n))
    print(time() - t)

    n = mpz(9832748274329847239**11)
    t = time()
    print(perfect_power(n))
    print(time() - t)