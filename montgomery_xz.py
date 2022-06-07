def xADD(P, Q, PQ, E):
    """
    Source: 1987 Montgomery 
    "Speeding the Pollard and elliptic curve methods of factorization", 
    
    4 Mult, 2 Square, 4 add 
    """
    N, _, A24 = E
    (PX, PZ), (QX, QZ), (PQX, PQZ) = P, Q, PQ

    A  = PX + PZ
    B  = PX - PZ
    C  = QX + QZ
    D  = QX - QZ
    DA = D * A
    CB = C * B
    X5 = (DA+CB)**2
    X5 = PQZ * X5
    Z5 = (DA-CB)**2
    Z5 = Z5 * PQX

    X5 = X5 % N
    Z5 = Z5 % N
    return (X5, Z5)

def xDBL(P, E):
    """
    1987 Montgomery 
    "Speeding the Pollard and elliptic curve methods of factorization",

    3 Mult, 2 Square, 4 Add 
    """
    N, _, A24 = E
    PX, PZ = P

    A  = PX + PZ
    AA = A**2
    B  = PX - PZ
    BB = B**2
    C  = AA - BB
    X3 = AA*BB
    Z3 = A24 * C
    Z3 = BB + Z3
    Z3 = C * Z3

    X3 = X3 % N
    Z3 = Z3 % N
    return (X3, Z3)

def xDBLADD(P, Q, PQ, E):
    """
    SIKE Specification Alg. 5
    https://sike.org/files/SIDH-spec.pdf

    7 Mult, 4 Sqaure, 8 Add
    """
    N, _, A24 = E
    (PX, PZ), (QX, QZ), (PQX, PQZ) = P, Q, PQ

    t0  = PX + PZ
    t1  = PX - PZ
    XPQ = QX + QZ
    t2  = QX - QZ
    X2P = t0 * t0
    t0  = t0 * t2
    Z2P = t1 * t1
    t1  = t1 * XPQ
    t2  = X2P - Z2P
    X2P = X2P * Z2P
    XPQ = A24 * t2
    ZPQ = t0 - t1
    Z2P = XPQ + Z2P
    XPQ = t0 + t1
    Z2P = Z2P * t2
    ZPQ = ZPQ * ZPQ
    XPQ = XPQ * XPQ
    ZPQ = PQX * ZPQ
    XPQ = PQZ * XPQ

    X2P = X2P % N
    Z2P = Z2P % N
    XPQ = XPQ % N
    ZPQ = ZPQ % N

    return (X2P, Z2P), (XPQ, ZPQ)

def xMUL(P, n, E):
    if n == 0:
        return (0, 0)
    R0, R1 = P, xDBL(P, E)
    for i in range(n.bit_length() - 2, -1, -1):
        if n & (1 << i) == 0:
            R0, R1 = xDBLADD(R0, R1, P, E)
        else:
            R1, R0 = xDBLADD(R0, R1, P, E)
    return R0