# gfedensity.sage
#
# Functions for computing the probability that Ax^l + By^m + Cz^n = 0
# has a primitive p-adic solution.
#
# See paper "On p-adic solubility of Ax^l + By^m + Cz^n = 0" by Christopher Keyes and Andrew Kobin
# ***URL***
#
# Authors: Christopher Keyes, Andrew Kobin
# Updated: 20 February 2026

# Set variables
var('p')             # p stands in for a prime
var('gcd_pmin1_l_m') # parameter for gcd(p-1,l,m), which will be substituted later
var('gcd_pmin1_l_n')
var('gcd_pmin1_m_n')

# solve_aux_x
#     Solves for the auxiliary probabilities rho_10c^(x), rho_1b0^(x), and analogous sigmas
#     See Lemma 4.5
# INPUTS:
#     * l,m,n: positive integer exponents
#     * gcd_pmin1_l_m, gcd_pmin1_l_n,gcd_pmin1_m_n: variables standing in for gcd(p-1,l,m), etc
#         (can also be a positive integer)
# OUTPUTS
#     * L_r_10c_vals: list of rho_10c^(x) values, as symbolic rational function expressions in
#         p (and gcd_pmin1_l_m, etc). Indexed by c-value, so L_r_10c_vals[c] = rho_10c^(x)
#     * L_r_1b0_vals, L_s_10c_vals, L_s_1b0_vals: analogous lists
def solve_aux_x(l,m,n,gcd_pmin1_l_m,gcd_pmin1_l_n,gcd_pmin1_m_n):
    # define relevant auxiliary probabilities
    var('r_100')
    var('s_100')

    # make lists of rho_10c^(x) for c < n
    L_r_10c = [r_100]
    L_s_10c = [s_100]
    for c in [1 .. n-1]:
        vstr = 'r_10'
        vstr = vstr + str(c)
        L_r_10c.append(var(vstr)) 
        vstr = 's_10'
        vstr = vstr + str(c)
        L_s_10c.append(var(vstr))    

    # make lists of rho_1b0^(x) for b < m
    L_r_1b0 = [r_100]
    L_s_1b0 = [s_100]
    for b in [1 .. m-1]:
        vstr = 'r_1'
        vstr = vstr + str(b)
        vstr = vstr + '0'
        L_r_1b0.append(var(vstr))  
        vstr = 's_1'
        vstr = vstr + str(b)
        vstr = vstr + '0'
        L_s_1b0.append(var(vstr))

    # build variable lists, initialize relation lists
    r_var_list = L_r_10c + [v for v in L_r_1b0 if v != r_100]
    s_var_list = L_s_10c + [v for v in L_s_1b0 if v != s_100]
    r_rel_list = []
    s_rel_list = []

    # relations for r_100, s_100
    if m == n:
        r_rel_list.append((r_100 == 1/gcd_pmin1_m_n + (1 - 1/gcd_pmin1_m_n)*((p-1)/p^m + 1/p^m*s_100)))    
        s_rel_list.append((s_100 == (p-1)/p^m + 1/p^m*s_100))
    elif m < n:
        r_rel_list.append((r_100 == 1/gcd_pmin1_m_n + (1 - 1/gcd_pmin1_m_n)*((p-1)/p^m/gcd_pmin1_l_m + 1/p^m*L_s_10c[n-m])))
        s_rel_list.append((s_100 == (p-1)/p^m/gcd_pmin1_l_m + 1/p^m*L_s_10c[n-m]))
    elif m > n:
        r_rel_list.append((r_100 == 1/gcd_pmin1_m_n + (1 - 1/gcd_pmin1_m_n)*((p-1)/p^n/gcd_pmin1_l_n + 1/p^n*L_s_1b0[m-n])))
        s_rel_list.append((s_100 == (p-1)/p^n/gcd_pmin1_l_n + 1/p^n*L_s_1b0[m-n]))

    # relations for r_10c and s_10c
    for c in [1 .. n-1]:
        if c == m:
            r_rel_list.append((L_r_10c[c] == (p-1)/p^m + 1/p^m*r_100))
            s_rel_list.append((L_s_10c[c] == (p-1)/p^m + 1/p^m*s_100))
        elif c < m:
            r_rel_list.append((L_r_10c[c] == (p-1)/p^c/gcd_pmin1_l_n + 1/p^c*L_r_1b0[m-c]))
            s_rel_list.append((L_s_10c[c] == (p-1)/p^c/gcd_pmin1_l_n + 1/p^c*L_s_1b0[m-c]))
        elif c > m:
            r_rel_list.append((L_r_10c[c] == (p-1)/p^m/gcd_pmin1_l_m + 1/p^m*L_r_10c[c-m]))
            s_rel_list.append((L_s_10c[c] == (p-1)/p^m/gcd_pmin1_l_m + 1/p^m*L_s_10c[c-m]))

    # relations for r_1b0 and s_1b0
    for b in [1 .. m-1]:
        if b == n:
            r_rel_list.append((L_r_1b0[b] == (p-1)/p^n + 1/p^n*r_100))
            s_rel_list.append((L_s_1b0[b] == (p-1)/p^n + 1/p^n*s_100))
        elif b < n:
            r_rel_list.append((L_r_1b0[b] == (p-1)/p^b/gcd_pmin1_l_m + 1/p^b*L_r_10c[n-b]))
            s_rel_list.append((L_s_1b0[b] == (p-1)/p^b/gcd_pmin1_l_m + 1/p^b*L_s_10c[n-b]))
        elif b > n:
            r_rel_list.append((L_r_1b0[b] == (p-1)/p^n/gcd_pmin1_l_n + 1/p^n*L_r_1b0[b-n]))
            s_rel_list.append((L_s_1b0[b] == (p-1)/p^n/gcd_pmin1_l_n + 1/p^n*L_s_1b0[b-n]))

    # solve system
    rel_list = r_rel_list + s_rel_list
    var_list = r_var_list + s_var_list
    S = solve(rel_list, var_list)
    
    # extract values from S
    L_r_10c_vals = [s.right() for s in S[0] if s.left() in L_r_10c]
    L_r_1b0_vals = [s.right() for s in S[0] if s.left() in L_r_1b0]
    L_s_10c_vals = [s.right() for s in S[0] if s.left() in L_s_10c]
    L_s_1b0_vals = [s.right() for s in S[0] if s.left() in L_s_1b0]
    
    return L_r_10c_vals, L_r_1b0_vals, L_s_10c_vals, L_s_1b0_vals

# get_rAB
#     Solves for rho_AB, given list of appropriate auxiliary probabilities
#     See Lemma 4.4
# INPUTS:
#     * L_r_10c_x: list of auxiliary probabilities rho_10c^(x)
#     * L_r_01c_y: list of auxiliary probabilities rho_01c^(y)
#     * n: positive integer exponent on z
#     * gcd_pmin1_l_m: variable standing in for gcd(p-1,l,m)
#         (can also be a positive integer)
# OUTPUTS
#     * rAB: symbolic rational function expression for rho_AB in variables
#         p and gcd_pmin1_l_m
def get_rAB(L_r_10c_x,L_r_01c_y,n,gcd_pmin1_l_m):
    # initialize with terms coming from min(a,b,n) = n
    rAB = (p-1)^2/p^(2*n) + (p-1)/p^(2*n)*(L_r_10c_x[0] + L_r_01c_y[0])
    
    # add terms coming from min(a,b,n) = i < n
    for i in [1 .. n-1]:
        r_10c_x = L_r_10c_x[n-i]
        r_01c_y = L_r_01c_y[n-i]
        rAB = rAB + (p-1)^2/p^(2*i)/gcd_pmin1_l_m + (p-1)/p^(2*i)*(r_10c_x + r_01c_y)

    return rAB/(1 - 1/p^(2*n))

# get_probs
#     Solves for rho, and other related probabilities for generic primes p,
#     i.e. those which do not divide gcd(l,m), gcd(l,n), gcd(m,n) for which 
#     rho_000^(x) = rho_000^(y) = rho_000^(z) = 1. See Theorem 5.2
# INPUTS:
#     * l,m,n: positive integer exponents
#     * (optional) return_all: true or false (default false - you probably don't need all)
# OUTPUTS
#     * r, rA, rB, rC, rAB, rAC, rBC (if return_all = false)
#     * also returns an excessive number of lists of auxiliary probabilties, if return_all = true
# NOTE: to be able to use get_rAB to calculate rAC and rBC, we need to run solve_aux_x 
#     for all permuatation of the variables/exponents. This is a little silly, but 
#     it works! We can write less code at the expense of having to keep track of
#     which gcd variable goes with which calculation.
def get_probs(l,m,n,return_all = false):
    # compute auxiliary probabilities for all permutations of variables
    L_r_10c_lmn, L_r_1b0_lmn, L_s_10c_lmn, L_s_1b0_lmn = solve_aux_x(l,m,n,gcd_pmin1_l_m,gcd_pmin1_l_n,gcd_pmin1_m_n)
    L_r_10c_lnm, L_r_1b0_lnm, L_s_10c_lnm, L_s_1b0_lnm = solve_aux_x(l,n,m,gcd_pmin1_l_n,gcd_pmin1_l_m,gcd_pmin1_m_n)
    L_r_10c_mln, L_r_1b0_mln, L_s_10c_mln, L_s_1b0_mln = solve_aux_x(m,l,n,gcd_pmin1_l_m,gcd_pmin1_m_n,gcd_pmin1_l_n)
    L_r_10c_mnl, L_r_1b0_mnl, L_s_10c_mnl, L_s_1b0_mnl = solve_aux_x(m,n,l,gcd_pmin1_m_n,gcd_pmin1_l_m,gcd_pmin1_l_n)
    L_r_10c_nlm, L_r_1b0_nlm, L_s_10c_nlm, L_s_1b0_nlm = solve_aux_x(n,l,m,gcd_pmin1_l_n,gcd_pmin1_m_n,gcd_pmin1_l_m)
    L_r_10c_nml, L_r_1b0_nml, L_s_10c_nml, L_s_1b0_nml = solve_aux_x(n,m,l,gcd_pmin1_m_n,gcd_pmin1_l_n,gcd_pmin1_l_m)

    # extract rho_A, rho_B, rho_C
    rA = L_r_10c_lmn[0]
    rB = L_r_10c_mln[0]
    rC = L_r_10c_nlm[0]
    
    # calculate rho_AB, rho_AC, rho_BC
    rAB = get_rAB(L_r_10c_lmn, L_r_10c_mln, n, gcd_pmin1_l_m)
    rAC = get_rAB(L_r_10c_lnm, L_r_10c_nlm, m, gcd_pmin1_l_n)
    rBC = get_rAB(L_r_10c_mnl, L_r_10c_nml, l, gcd_pmin1_m_n)

    # calculate rho
    r = ((p-1)^3/p^3 + (p-1)^2/p^3*(rA + rB + rC) + (p-1)/p^3*(rAB + rAC + rBC))/(1-1/p^3)
    
    # return literally everything (you probably don't need this)
    if return_all:
        return r, rA, rB, rC, rAB, rAC, rBC, L_r_10c_lmn, L_r_1b0_lmn, L_s_10c_lmn, L_s_1b0_lmn, L_r_10c_mln, L_r_1b0_mln, L_s_10c_mln, L_s_1b0_mln, L_r_10c_nlm, L_r_1b0_nlm, L_s_10c_nlm, L_s_1b0_nlm
    
    return r, rA, rB, rC, rAB, rAC, rBC

# validity_range
#     returns the last prime (not dividing l,m,n) for which the rational function expressions
#     does NOT (necessarily) compute the p-adic probability. This is computed via Theorem 1.1
#     (and Lemmas 2.8, 3.6), and can be improved in many situation.
# INPUTS:
#     * l,m,n: positive integer exponents
# OUTPUTS
#     * last_prime: the last prime (not dividing l,m,n) for which the result is not guaranteed
def validity_range(l,m,n):
    g = gcd(gcd(l,m),n)
    
    if g == 1:
        return 2
    elif g == 2:
        return 5
    
    stop_search = false
    n = 1
    bd_old = n + 1 - (g-2)*(g-1)/2*floor(2*sqrt(n)) - 3*g
    last_prime = 1
    
    # increment n until this bound drops, but is still positive
    while not stop_search:
        n = n+1
        bd_new = n + 1 - (g-2)*(g-1)/2*floor(2*sqrt(n)) - 3*g
        # printing for debugging
#         print(n,bd_new)
        
        if bd_new > 0:
            if bd_new <= bd_old:
                stop_search = true
        
        # bound is nonpositive
        else:
            # if n prime, record actual bound (incorporating gcd) and update last_prime if nonpositive
            if n.is_prime():
                bd_actual = n+1 - (g-2)*(g-1)/2*floor(2*sqrt(n)) - 3*gcd(n-1,g)
                if bd_actual <= 0:
                    last_prime = n
                
        bd_old = bd_new
        
    return last_prime

# validity_range_diagonal
#     returns the last prime (not dividing n) for which the rational function expressions
#     do NOT (necessarily) compute the p-adic probability when (l,m,n) = (n,n,n). This 
#     takes advantage of improvements to Theorem 1.1 in the diagonal case.
# INPUTS:
#     * n: positive integer exponent
# OUTPUTS
#     * last_prime: the last prime (not dividing n) for which the result is not guaranteed
def validity_range_diagonal(n):
    # if p is big enough, Hirakawa--Kanamura (remark 2.4) show validity 
    plist = [p for p in [2 .. (n-1)^2*(n-2)^2] if p.is_prime() and mod(n,p) != 0]    
    
    # initialize
    last_prime = 1
    
    for p in plist:
        # if p is big enough (Hasse--Weil bound) or gcd(p-1,n) = 1, we are guaranteed points
        if p+1 - (n-1)*(n-2)/2*floor(2*sqrt(p)) <= 0 and gcd(p-1,n) != 1:
            last_prime = p
            
    return last_prime

# ptlists
#     returns list of nontrivial (x,y,z) solving Ax^l + By^m + Cz^n = 0 over finite field Fp
#     and lists of those with x,y,z nonzero
# INPUTS:
#     * l,m,n: positive integer exponents
#     * A,B,C: elements of finite field Fp (integers will also work)
#     * p: a prime
# OUTPUTS
#     * ptlist: list of (x,y,z) satisfying the equation
#     * ptlist_x, ptlist_y, ptlist_z: sublists with x (resp. y, z) nonzero
def ptlists(l,m,n,A,B,C,p):
    F = FiniteField(p)
    xyzlist = [(x,y,z) for (x,y,z) in tuples(F,3) if x.is_unit() or y.is_unit() or z.is_unit()]
    ptlist = [(x,y,z) for (x,y,z) in xyzlist if A*x^l + B*y^m + C*z^n == 0]
    ptlist_x = [(x,y,z) for (x,y,z) in ptlist if x.is_unit()]
    ptlist_y = [(x,y,z) for (x,y,z) in ptlist if y.is_unit()]
    ptlist_z = [(x,y,z) for (x,y,z) in ptlist if z.is_unit()]
    return ptlist, ptlist_x, ptlist_y, ptlist_z

# rho0_xyz
#     computes rho0 (and rho0^(x), rho0^(y), rho0^(z)) for primes p dividing no two of l,m,n
# INPUTS:
#     * l,m,n: positive integer exponents
#     * p: a prime which we assume divides no two of l,m,n (this is asserted)
# OUTPUTS:
#     * rho0: the ratio of A,B,C units in Fp with solutions to the total
#     * rho0^(x), rho0^(y), rho0^(z): the analogous ratios
# NOTE: 
#     * this ratio computes rho0 because so long as p doesn't divide any two of l,m,n, since 
#         any solution (x,y,z) has p dividing at most one coordinate, we must be able to lift via
#         Hensel's lemma in one coordinate to produce a Zp-solution. Conversely, if no Fp-solution
#         exists, no Zp-solution does.
#     * this is not optimized at all. In fact it is quite inefficient, as it computes the list
#         of ALL Fp-solutions to the equation, rather than just finding one and quitting. 
#         As such, NOT RECOMMENDED for p much bigger than very tiny. 
def rho0_xyz(l,m,n,p):
    assert mod(gcd(l,m),p) != 0 and mod(gcd(l,n),p) != 0 and mod(gcd(m,n),p) != 0
    
    F = FiniteField(p)
    Fnz = [a for a in F if a.is_unit()] # nonzero elements of F
    
    # initialize counts
    ct = 0
    ctpt = 0
    ctx = 0
    cty = 0
    ctz = 0
    
    # suffices to scale so A=1 and count over B,C units in F
    for B,C in tuples(Fnz,2):
        ct = ct + 1
        
        # overkill, but we have a routine to get lists of points
        ptlist, ptlist_x, ptlist_y, ptlist_z = ptlists(l,m,n,1,B,C,p)
        if len(ptlist) > 0:
            ctpt = ctpt + 1
        if len(ptlist_x) > 0:
            ctx = ctx + 1
        if len(ptlist_y) > 0:
            cty = cty + 1
        if len(ptlist_z) > 0:
            ctz = ctz + 1
            
    return ctpt/ct, ctx/ct, cty/ct, ctz/ct        

# valid_small_p
#     checks if conditions of Theorem 5.2 are valid for primes smaller than the bounds given by 
#     validity_range (or other methods). NOT RECOMMENDED for primes p bigger than very tiny
# INPUTS:
#     * l,m,n: positive integer exponents
#     * last_prime: a threshhold for which we know the formulae are valid for p > last_prime
#         (e.g. that produced by validity_range)
# OUTPUTS:
#     * valid_plist: list of primes between 2 and last_prime for which p doesn't divide any two
#         of l,m,n and rho0^(x) = rho0^(y) = rho0^(z) = 1
# NOTE: unoptimized and inefficient - see comments for rho0_xyz
def valid_small_p(l,m,n,last_prime):
    plist = [p for p in [2 .. last_prime] if p.is_prime()]
    valid_plist = []
    
    for p in plist:
        # check if p divides any two of l,m,n
        if not((mod(gcd(l,m),p) == 0 and mod(gcd(l,n),p) == 0) or (mod(gcd(l,m),p) == 0 and mod(gcd(m,n),p) == 0) or (mod(gcd(l,n),p) == 0 and mod(gcd(m,n),p) == 0)):
            # compute rho0's
            r,rx,ry,rz = rho0_xyz(l,m,n,p)
            # if all 1, then it's valid!
            if r == 1 and rx == 1 and ry == 1 and rz == 1:
                valid_plist.append(p)                
    
    return valid_plist

# possible_gcds
#     returns possible values of i=gcd(p-1,l,m), j=gcd(p-1,l,n), k=gcd(p-1,m,n)
# INPUTS:
#     * l,m,n: positive integer exponents
# OUTPUTS:
#     * L_ijk: list of possible triples (i,j,k)
def possible_gcds(l,m,n):
    glm = gcd(l,m)
    gln = gcd(l,n)
    gmn = gcd(m,n)
    
    L_ijk = []
    for i in divisors(glm):
        for j in divisors(gln):
            for k in divisors(gmn):
                if gcd(i,gln) in divisors(j) and gcd(i,gmn) in divisors(k) and gcd(j,glm) in divisors(i) and gcd(j,gmn) in divisors(k) and gcd(k,glm) in divisors(i) and gcd(k,gln) in divisors(j):
                    L_ijk.append((i,j,k))
    
    return L_ijk

# print_ratfunc
#     prints a rational function with specialization i=gcd(p-1,l,m), j=gcd(p-1,l,n), k=gcd(p-1,m,n)
# INPUTS:
#     * R: a rational function in variables p and gcd_pmin1_l_m, gcd_pmin1_l_n, gcd_pmin1_m_n
#     * i: positive integer value for gcd(p-1,l,m)
#     * j: positive integer value for gcd(p-1,l,n)
#     * k: positive integer value for gcd(p-1,m,n)
#     * texprint (optional): true or false (prints latex if true, numerator/denominator if false)
# NOTE: it's nice to have some flexibility with R here, so we can print things like 1-rho, etc.
def print_ratfunc(R,i,j,k,texprint=false):
    
    if texprint:
        print("\\left(")
        print(latex(R(gcd_pmin1_l_m = i, gcd_pmin1_l_n = j, gcd_pmin1_m_n = k).expand().numerator()))
        print("\\right)\\Big/\\left(")
        print(latex(R(gcd_pmin1_l_m = i, gcd_pmin1_l_n = j, gcd_pmin1_m_n = k).expand().denominator()))
        print("\\right)")
    else:    
        print(R(gcd_pmin1_l_m = i, gcd_pmin1_l_n = j, gcd_pmin1_m_n = k).expand().numerator().factor())
        print(R(gcd_pmin1_l_m = i, gcd_pmin1_l_n = j, gcd_pmin1_m_n = k).expand().denominator().factor())
        
# print_all_ratfunc
#     prints all specializations i=gcd(p-1,l,m), j=gcd(p-1,l,n), k=gcd(p-1,m,n) for some rational
#     function
# INPUTS:
#     * l,m,n: positive integer exponents
#     * R: a rational function in variables p and gcd_pmin1_l_m, gcd_pmin1_l_n, gcd_pmin1_m_n
#     * Rname: a string name for the rational function to be printed
#     * texprint (optional): true or false (prints latex if true, numerator/denominator if false)
# NOTE: it's nice to have some flexibility with R here, so we can print things like 1-rho, etc.
def print_all_ratfunc(l,m,n,R,Rname="no name provided",texprint=false):
    # get list of i,j,k values
    L_ijk = possible_gcds(l,m,n)
    
    print("l,m,n =", l,m,n)
    print("printing ", Rname)
    # optional: printing validity range
#     last_prime = validity_range(l,m,n)
#     print("p >", last_prime, "and p doesn't divide lmn")
#     print("also valid for p in", valid_small_p(l,m,n,last_prime))
    print("possible i,j,k =", L_ijk, "\n")
    
    # print all rational function specializations
    for ijk in L_ijk:
        i,j,k = ijk
        print("i,j,k =", i,j,k)
        print_ratfunc(R,i,j,k,texprint=texprint)
        print("")
