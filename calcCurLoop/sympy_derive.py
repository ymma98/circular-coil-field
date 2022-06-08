
"""
ymma, ymma98@qq.com
2022.06.08
"""
import sympy as sym

"""
ref:
[1] Jackson, Thrid edition, eq. (5.38)
[2] Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop
    James Simpson et al
notice: In Jackson's book, eq.(5.37) is wrong. For elliptic function on the numerator,
        they should be E(k**2) and K(k**2), but not the E(k) and K(k)
"""

r, th, mu, I, a = sym.symbols('r, theta, mu, I, a')
k = sym.sqrt((4*a*r*sym.sin(th)) / (a**2 + r**2 + 2*a*r*sym.sin(th)))
aphi = mu / (sym.pi * 4) * (4*I*a) / (sym.sqrt(a**2 + r**2 + 2*a*r*sym.sin(th))) * ((2-k**2)*sym.elliptic_k(k**2) - 2*sym.elliptic_e(k**2)) / k**2
Br = 1/(r*sym.sin(th)) * sym.diff(sym.sin(th)*aphi, th)
Bth = -1/r * sym.diff(r*aphi, r)


print("*****"*2 + "B_r" + "*****"*2)
print("*****"*2)
print("*****" + "pprint version" + "*****")
sym.pprint(sym.simplify(Br))
print("*****" + "latex version" + "*****")
print(sym.latex(sym.simplify(Br)))
print("*****" + "python version" + "*****")
print(sym.python(sym.simplify(Br)))
print("*****"*2)

print('\n\n\n')
print("*****"*2 + "B_theta" + "*****"*2)
print("*****"*2)
print("*****" + "pprint version" + "*****")
sym.pprint(sym.simplify(Bth))
print("*****" + "latex version" + "*****")
print(sym.latex(sym.simplify(Bth)))
print("*****" + "python version" + "*****")
print(sym.python(sym.simplify(Bth)))
print("*****"*2)
