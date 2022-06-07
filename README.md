# Lenstra ECM

Implmentation of Lenstra's ECM factorisation algorithm following Pomerance's "Prime Numbers: A Computational Perspective", Algorithm 7.4.4 (Inversionless ECM). 

## ECM Factoring

:construction: Under construction :construction:

Maybe I could use the nice Latex support GitHub now has to give an overview of how this algorithm works... :eyes:

## cProfile

Running `cProfile` very clearly shows that the bulk of the computation happens in performing `xDBLADD`.


```py
         807298 function calls in 3.361 seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    3.361    3.361 <string>:1(<module>)
        1    0.002    0.002    3.361    3.361 factor.py:5(factor)
       35    0.000    0.000    0.001    0.000 lenstra_ecm.py:30(find_curve)
        8    1.608    0.201    3.345    0.418 lenstra_ecm.py:52(lenstra_ecm)
    43876    0.087    0.000    0.087    0.000 montgomery_xz.py:1(xADD)
    46837    0.063    0.000    0.063    0.000 montgomery_xz.py:26(xDBL)
   534198    1.358    0.000    1.358    0.000 montgomery_xz.py:50(xDBLADD)
    46781    0.209    0.000    1.633    0.000 montgomery_xz.py:87(xMUL)
       35    0.000    0.000    0.000    0.000 random.py:239(_randbelow_with_getrandbits)
       35    0.000    0.000    0.000    0.000 random.py:292(randrange)
       35    0.000    0.000    0.000    0.000 random.py:366(randint)
        1    0.013    0.013    0.014    0.014 trial_division.py:4(trial_division)
     4160    0.003    0.000    0.003    0.000 {built-in method _bisect.bisect_left}
     4168    0.003    0.000    0.003    0.000 {built-in method _bisect.bisect_right}
      105    0.000    0.000    0.000    0.000 {built-in method _operator.index}
        1    0.000    0.000    3.361    3.361 {built-in method builtins.exec}
        1    0.000    0.000    0.000    0.000 {built-in method builtins.min}
       70    0.000    0.000    0.000    0.000 {built-in method builtins.pow}
        2    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.ceil}
       63    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.gcd}
        9    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.is_prime}
        2    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.isqrt}
        8    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.sqrt}
    10680    0.002    0.000    0.002    0.000 {built-in method math.log}
        9    0.005    0.001    0.005    0.001 {built-in method primesieve._primesieve.primes}
    69312    0.005    0.000    0.005    0.000 {method 'append' of 'list' objects}
    46816    0.003    0.000    0.003    0.000 {method 'bit_length' of 'int' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
       48    0.000    0.000    0.000    0.000 {method 'getrandbits' of '_random.Random' objects}
```

Thoughts:

- Would implementing Montgomery Arithmetic for `xDBLADD` speed it up, or slow it down due to the overhead? 11 Multiplications isn't so much...
- Zimmerman has a very different approach for Stage Two, which *also* uses affine coordinates for Stage Two. I need to learn how this works, implement and then compare. Currently Stage One is much shorter than Stage Two, and from the literature I think bounds should be balanced so they take approximately the same time?
- Would be nice to include a $P-1$ factoring if ECM fails, just incase we get lucky
- I currently don't have any check that $N \neq x^k$. Perfect powers could be tested after trial division and before ECM?
- Rather than implementing all of this in Rust/C, would there be too much FFI overhead to shift the group operations into C and then hook these back into Python?

