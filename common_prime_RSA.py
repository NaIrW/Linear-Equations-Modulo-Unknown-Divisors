import itertools
from subprocess import check_output
from re import findall


def flatter(M):
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(rb"-?\d+", ret)))


def solve(N, e, t=5):
    gamma = 0.42
    beta = 0.25

    PR.<x, y> = ZZ[]
    E = inverse_mod(e, (N - 1) // 2)
    f1 = E - x
    f2 = N - y

    bounds = (2^250, 2^500)
    U = N - 1
    
    polynomials = Sequence([], PR)
    for i1, i2 in itertools.product(range(ceil(gamma * t / beta)), range(ceil(gamma * t / 0.5))):
        if beta * i1 + 0.5 * i2 <= gamma * t:
            polynomials.append(f1^i1 * f2^i2 * U^max(0, t - i1 - 2 * i2))
    polynomials.sort()
    B, monomials = polynomials.coefficient_matrix()
    monomials = vector(monomials)

    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)

    print(f'dim(L): {len(monomials)}')
    B = flatter(B)
    print('flatter done')

    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1 / factor)
    B *= monomials

    y = y.change_ring(QQ)
    PR.<q> = ZZ[]
    for pol1, pol2 in itertools.combinations(B, 2):
        rr = pol1.resultant(pol2, y)(q, q)
        solx = rr.roots()
        if solx:
            return solx[0][0]
                

N = 0xbe9ccc83003bedf45421b58377b946f87dfd85be82124dc5d732070d77ef68e0231c3f34dc803a8984de0573db6d83ccea0bd53a885059a10cfa3764c658c4d42c5fa90ecad8573fff8f2c41e513278c59121e42ad83310fb22b4d20e7ada42c76f08891f38c92a1b1aac712bfa7d717a4c4802ed023f12c768972ca1b
e = 0x5dc97ed7250e57ce6fac4f57885c0538b1ea540fbaca79730470b6b990f7e861adc4c5fee3acdcd9ae9a2834b606ddfae01ade33edfa96a47a0ffc0036a4497a84c38b7cdac20c38f
d = solve(N, e, t=10)
print(f'd = {hex(d)}')
print('DubheCTF{{{:x}}}'.format(ZZ(d) & ZZ(2^128 - 1)))

"""
d = 0x2e2c8c829244854b305b75cce11f62eb896a5fef7abec06cd2e6256be4ba40b
DubheCTF{b896a5fef7abec06cd2e6256be4ba40b}
"""
