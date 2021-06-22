# -*- coding: utf-8

# License
# Copyright

r"""
Diagonals in the "double-(2g+1)-gon staircase"

EXAMPLES::

    sage: # load('hecke_orbits.sage')
    sage: basic_diagonals(1)
    [(1, 0), (1, 1)]
    sage: basic_diagonals(2)
    [(1, 0), (a, 1), (a, a)]
    sage: basic_diagonals(3)
    [(1, 0), (a, 1), (a^2 - 1, a), (a^2 - 1, a^2 - 1)]
    sage: basic_diagonals(4)
    [(1, 0), (a, 1), (a^2 - 1, a), (a + 1, a^2 - 1), (a + 1, a + 1)]
    sage: basic_diagonals(5)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^4 - 3*a^2 + 1, a^4 - 3*a^2 + 1)]
    sage: basic_diagonals(6)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^5 - 4*a^3 + 3*a, a^4 - 3*a^2 + 1),
     (a^5 - 4*a^3 + 3*a, a^5 - 4*a^3 + 3*a)]
    sage: basic_diagonals(7)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (-a^3 + a^2 + 4*a, a^3 - 2*a),
     (a^3 - 2*a + 1, -a^3 + a^2 + 4*a),
     (a^2 + a - 1, a^3 - 2*a + 1),
     (a^2 + a - 1, a^2 + a - 1)]
    sage: basic_diagonals(8)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^5 - 4*a^3 + 3*a, a^4 - 3*a^2 + 1),
     (a^6 - 5*a^4 + 6*a^2 - 1, a^5 - 4*a^3 + 3*a),
     (a^7 - 6*a^5 + 10*a^3 - 4*a, a^6 - 5*a^4 + 6*a^2 - 1),
     (a^7 - 6*a^5 + 10*a^3 - 4*a, a^7 - 6*a^5 + 10*a^3 - 4*a)]
    sage: basic_diagonals(9)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^5 - 4*a^3 + 3*a, a^4 - 3*a^2 + 1),
     (a^6 - 5*a^4 + 6*a^2 - 1, a^5 - 4*a^3 + 3*a),
     (a^7 - 6*a^5 + 10*a^3 - 4*a, a^6 - 5*a^4 + 6*a^2 - 1),
     (a^8 - 7*a^6 + 15*a^4 - 10*a^2 + 1, a^7 - 6*a^5 + 10*a^3 - 4*a),
     (a^8 - 7*a^6 + 15*a^4 - 10*a^2 + 1, a^8 - 7*a^6 + 15*a^4 - 10*a^2 + 1)]
    sage: basic_diagonals(10)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^5 - 4*a^3 + 3*a, a^4 - 3*a^2 + 1),
     (-a^5 + a^4 + 6*a^3 - 2*a^2 - 8*a - 2, a^5 - 4*a^3 + 3*a),
     (a^5 - 4*a^3 + 3*a + 1, -a^5 + a^4 + 6*a^3 - 2*a^2 - 8*a - 2),
     (a^4 - 3*a^2 + a + 1, a^5 - 4*a^3 + 3*a + 1),
     (a^3 + a^2 - 2*a - 1, a^4 - 3*a^2 + a + 1),
     (a^3 + a^2 - 2*a - 1, a^3 + a^2 - 2*a - 1)]
    sage: basic_diagonals(11)
    [(1, 0),
     (a, 1),
     (a^2 - 1, a),
     (a^3 - 2*a, a^2 - 1),
     (a^4 - 3*a^2 + 1, a^3 - 2*a),
     (a^5 - 4*a^3 + 3*a, a^4 - 3*a^2 + 1),
     (a^6 - 5*a^4 + 6*a^2 - 1, a^5 - 4*a^3 + 3*a),
     (a^7 - 6*a^5 + 10*a^3 - 4*a, a^6 - 5*a^4 + 6*a^2 - 1),
     (a^8 - 7*a^6 + 15*a^4 - 10*a^2 + 1, a^7 - 6*a^5 + 10*a^3 - 4*a),
     (a^9 - 8*a^7 + 21*a^5 - 20*a^3 + 5*a, a^8 - 7*a^6 + 15*a^4 - 10*a^2 + 1),
     (a^10 - 9*a^8 + 28*a^6 - 35*a^4 + 15*a^2 - 1, a^9 - 8*a^7 + 21*a^5 - 20*a^3 + 5*a),
     (a^10 - 9*a^8 + 28*a^6 - 35*a^4 + 15*a^2 - 1, a^10 - 9*a^8 + 28*a^6 - 35*a^4 + 15*a^2 - 1)]

    sage: # load('hecke_orbits.sage')
    sage: d = long_diagonals(2, 9, 9)  # long time
    sage: pic = nf_point2d(d) + line2d(((5,0),(5,5),(0,5)))  # long time
    sage: ppic = pic + sum(line2d(((0,0),p)) for p in d)  # long time
    sage: ppic.show(aspect_ratio=1,figsize=(20,20))  # long time

    sage: # load('hecke_orbits.sage')
    sage: xmax, ymax = 100, 100  # long time
    sage: xmax, ymax = 500, 500  # long time
    sage: d = long_diagonals(2, xmax, ymax, clip=True)  # long time
    sage: pic = nf_point2d(d) + line2d(((xmax, 0),(xmax, ymax), (0, ymax)))  # long time
    sage: pic.show(aspect_ratio=1,figsize=(20, 20))  # long time
    sage: nf_point2d(d).show(aspect_ratio=1, axes=False, figsize=(100,100))  # long time

    sage: # load('hecke_orbits.sage')
    sage: xmax, ymax = 100, 100  # long time
    sage: xmax, ymax = 500, 500  # long time
    sage: d2 = long_diagonals(2, xmax, ymax, clip=True)  # long time
    sage: nf_point2d(d2).show(aspect_ratio=1, axes=False, figsize=(xmax/20, ymax/20))  # long time
    sage: d3 = long_diagonals(3, xmax, ymax, clip=True)  # long time
    sage: nf_point2d(d3).show(aspect_ratio=1, axes=False, figsize=(xmax/20, ymax/20))  # long time
    sage: d4 = long_diagonals(4, xmax, ymax, clip=True)  # long time
    sage: nf_point2d(d4).show(aspect_ratio=1, axes=False, figsize=(xmax/20, ymax/20))  # long time
    sage: d5 = long_diagonals(5, xmax, ymax, clip=True)  # long time
    sage: nf_point2d(d5).show(aspect_ratio=1, axes=False, figsize=(xmax/20, ymax/20))  # long time
    sage: d6 = long_diagonals(6, xmax, ymax, clip=True)  # long time
    sage: nf_point2d(d6).show(aspect_ratio=1, axes=False, figsize=(xmax/20, ymax/20))  # long time

    sage: # load('hecke_orbits.sage')
    sage: gmin, gmax = 2, 8  # long time
    sage: xymax = 250  # long time
    sage: verbose = True  # long time
    sage: for g in sxrange(gmin, gmax+1):  # long time
    ....:     if verbose:
    ....:         print('g = {}'.format(g))
    ....:     p = nf_point2d(long_diagonals(g, xymax, xymax),
    ....:         aspect_ratio=1, figsize=(10, 10))
    ....:     p.save('long_diagonals_g_{}_xy_{}.svg'.format(g, xymax))
    ....:     p.save('long_diagonals_g_{}_xy_{}.png'.format(g, xymax))

    sage: d = long_diagonals(g=1, xmax=250, ymax=250)  # long time
    sage: dd = [vector(ZZ, v) for v in d if not (ZZ(v[0]) % 2) or not (ZZ(v[1]) % 2)]  # long time
    sage: p = nf_point2d(dd, aspect_ratio=1, figsize=(10, 10))  # long time
    sage: p.save('long_diagonals_g_{}_xy_{}.png'.format('oo', 250))  # long time
    sage: p.save('long_diagonals_g_{}_xy_{}.svg'.format('oo', 250))  # long time

.. NOTES:

    - The time-consuming steps are (1) plotting, (2) enumerating;
      not much point in speeding up the preliminary steps;
      and not clear whether this would speed them up:

          def xxx(g, n=None, a=None, K=None, V=None):
              if n is None or a is None or K is None or V is None:
                   g, n, a, K, V = structure(g)

.. TODO:

    - Read Koseleff-Rouillier-Tran,
      "On the Sign of a Trigonometric Expression"

    - Animate the staircase surfaces Stairs(g) as g grows::

          sage: gmax = 12
          sage: gmax = 40
          sage: animate(Stairs(g).show(aspect_ratio=1, axes=False)
          ....:         for g in range(3, gmax+1))  # long time

    - Animate "long diagonals of Stairs(g)" as g grows

          sage: gmax = 12
          sage: gmax = 40
          sage: animate(long_diagonals(g, 200, 200).show(
          ....:          aspect_ratio=1, figsize=(10, 10))
          ....:          for g in range(3, gmax+1))  # long time



          sage: longdiags = []
          sage: gmax = 27
          sage: xmax = 100
          sage: for g in sxrange(3, gmax+1):  # long time
          ....:     longdiags.append(nf_point2d(long_diagonals(g, xmax, xmax),
          ....:             aspect_ratio=1, figsize=(10, 10)))
          sage: animate(longdiags)  # long time

    - Animate "all diagonals of Stairs(g)" as g grows

          sage: gmax = 12
          sage: gmax = 40
          sage: longdiags = []  # long time
          sage: for g in sxrange(3, gmax+1):  # long time
          ....:     print(g, end=' ')
          ....:     longdiags.append(nf_point2d(long_diagonals(g, 30, 30),
          ....:                      aspect_ratio=1, figsize=(10, 10))
          sage: animate(longdiags, aspect_ratio=1, figsize=(10, 10)) for g in range(3,gmax+1)))

    - line2d of long diagonals

          sage: line2d(long_diagonals(3, 10, 10)).show(
          ....:          aspect_ratio=1, axes=False, figsize=(10, 10))  # long time

Finding holes.

    sage: def boxes(d, box_size=1, return_box=False):  # long time
    ....:     b = box_size
    ....:     from collections import defaultdict
    ....:     box = defaultdict(list)
    ....:     xmax = ymax = 2^4
    ....:     for x, y in d:
    ....:         box[(b*floor(x/b), b*floor(y/b))].append((x,y))
    ....:     return box
    sage: xmax = ymax = 2^10  # long time
    sage: d2 = long_diagonals(2, xmax, ymax, clip=True)  # long time
    sage: len(d2)  # long time
    572605
    sage: B = boxes(d2, box_size=4, return_box=True)  # long time
    sage: G = Graphics()  # long time
    sage: for b, l in box.iteritems():  # long time
    ....:     G += text('{0:02d}'.format(len(l)), b)
    sage: G.show(aspect_ratio=1, aspect_ratio=1, figsize=64)  # long time
"""
# from __future__ import print_function
from collections import defaultdict
from itertools import chain

from sage.modules.free_module import VectorSpace
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.qqbar import AA
from sage.symbolic.ring import SR
from sage.misc.cachefunc import cached_function

from rnb.cos_minpoly import *
from rnb.timings import *

@cached_function
def structure(g):
    """
    Return g, n, a, K, V

    OUTPUT:

    - ``g`` -- the genus
    - ``n`` (= 2 * g + 1) -- so that double-n-gon is related
    - ``a`` -- 2 * cos(pi/n), as a number field element
    - ``K`` -- the number field QQ(a)
    - ``V`` -- the vector field K^2

    .. TODO:

    - speed things up using Koseleff-Rouillier-Tran's paper
      "On the Sign of a Trigonometric Expression"

    EXAMPLES::

        sage: structure(2)
        (2, 5, a,
         Number Field in a with defining polynomial x^2 - x - 1,
         Vector space of dimension 2 over Number Field in a
         with defining polynomial x^2 - x - 1)
        sage: structure(3)
        (3, 7, a,
         Number Field in a with defining polynomial x^3 - x^2 - 2*x + 1,
         Vector space of dimension 2 over Number Field in a
         with defining polynomial x^3 - x^2 - 2*x + 1)
    """
    g = ZZ(g)
    n = 2 * g + 1
    # aa = AA((SR.pi()/n).cos()*2)
    # p = aa.minpoly()
    # K = NumberField(p, 'a', embedding=RDF(aa))
    aa = (RDF.pi()/n).cos()*2
    p = cos_minpoly(n, polygen(ZZ))
    K = NumberField(p, 'a', embedding=aa)
    a = K.gen()
    V = VectorSpace(K, 2)
    return g, n, a, K, V


@cached_function
def gnapKOVM(g):
    """
    Return g, n, a, p, K, O, V, M

    OUTPUT:

    - ``g`` -- the genus

    - ``n`` (= 2 * g + 1) -- so that double-n-gon is related

    - ``a`` -- ``2 * cos(pi/n)``, as a number field element

    - ``p`` -- the minimal polynomial of `a`

    - ``K`` -- the number field QQ(a)

    - ``O`` -- the ring of integers of `K`

    - ``V`` -- the vector space `K^2`

    - ``M`` -- the matrix space of 2 by 2 matrices over `K`

    .. TODO:

    - speed things up using Koseleff-Rouillier-Tran's paper
      "On the Sign of a Trigonometric Expression"

    EXAMPLES::

        sage: gnapKOVM(2)
        (2, 5, a,
         x^2 - x - 1,
         Number Field in a with defining polynomial x^2 - x - 1,
         Maximal order in Number Field in a with defining polynomial x^2 - x - 1,
         Vector space of dimension 2 over Number Field in a
         with defining polynomial x^2 - x - 1)
         Matrix space of 2 by 2 matrices over Number Field in a
         with defining polynomial x^2 - x - 1)
        (2, 5, a,
         x^2 - x - 1,
         Number Field in a with defining polynomial x^2 - x - 1
         with a = 1.618033988749895?,
         Maximal Order in Number Field in a with defining polynomial x^2 - x - 1
         with a = 1.618033988749895?,
         Vector space of dimension 2 over Number Field in a with defining
         polynomial x^2 - x - 1 with a = 1.618033988749895?,
         Full MatrixSpace of 2 by 2 dense matrices over Number Field in a
         with defining polynomial x^2 - x - 1 with a = 1.618033988749895?)
        sage: gnapKOVM(3)
        (3, 7, a,
         x^3 - x^2 - 2*x + 1,
         Number Field in a with defining polynomial x^3 - x^2 - 2*x + 1
         with a = 1.801937735804839?,
         Maximal Order in Number Field in a with defining polynomial x^3 - x^2 - 2*x + 1
         with a = 1.801937735804839?,
         Vector space of dimension 2 over Number Field in a with defining
         polynomial x^3 - x^2 - 2*x + 1 with a = 1.801937735804839?,
         Full MatrixSpace of 2 by 2 dense matrices over Number Field in a
         with defining polynomial x^3 - x^2 - 2*x + 1 with a = 1.801937735804839?)
    """
    g = ZZ(g)
    n = 2 * g + 1
    # aa = AA((SR.pi()/n).cos()*2)
    # p = aa.minpoly()
    # K = NumberField(p, 'a', embedding=RDF(aa))
    aa = (RDF.pi()/n).cos()*2
    p = cos_minpoly(n, polygen(ZZ))
    K = NumberField(p, 'a', embedding=aa)
    a = K.gen()
    O = K.ring_of_integers()
    V = VectorSpace(K, 2)
    M = MatrixSpace(K, 2)
    return g, n, a, p, K, O, V, M


@cached_function
def number_sequence(g, all_the_way=False):
    """
    Return number sequence giving coordinates of basic diagonals

    EXAMPLES::

        sage: number_sequence(2)
        (0, 1, a)
        sage: [number_sequence(g) for g in range(1, 6)]
        [(0, 1),
         (0, 1, a),
         (0, 1, a, a^2 - 1),
         (0, 1, a, a^2 - 1, a + 1),
         (0, 1, a, a^2 - 1, a^3 - 2*a, a^4 - 3*a^2 + 1)]
    """
    g, n, a, K, V = structure(g)
    s = [K.zero(), K.one()]
    while len(s) <= g:
        s.append(a*s[-1] - s[-2])
    if all_the_way:
        s.extend(s[::-1])
    return tuple(s)


@cached_function
def basic_diagonals(g, all_the_way=False):
    """
    Return the basic diagonals for the (odd n) "double-n-gon staircase".

    INPUT:

    - ``g`` -- integer: the genus

    - ``all_the_way`` -- boolean: whether to include basic diagonals
      of slope > 1.

    EXAMPLES::

        sage: basic_diagonals(2)
        [(1, 0), (a, 1), (a, a)]
        sage: [basic_diagonals(g) for g in range(1, 5)]
        [[(1, 0), (1, 1)],
         [(1, 0), (a, 1), (a, a)],
         [(1, 0), (a, 1), (a^2 - 1, a), (a^2 - 1, a^2 - 1)],
         [(1, 0), (a, 1), (a^2 - 1, a), (a + 1, a^2 - 1), (a + 1, a + 1)]]
        sage: [basic_diagonals(g, all_the_way=True) for g in range(1, 4)]
        [[(1, 0), (1, 1)],
         [(1, 0), (a, 1), (a, a)],
         [(1, 0), (a, 1), (a^2 - 1, a), (a^2 - 1, a^2 - 1)],
         [(1, 0), (a, 1), (a^2 - 1, a), (a + 1, a^2 - 1), (a + 1, a + 1)]]
    """
    g, n, a, K, V = structure(g)
    v = [V((K.one(), K.zero())), V((a, K.one()))]
    while len(v) <= g:
        v.append(a*v[-1] - v[-2])
    if all_the_way:
        v.extend(u[::-1] for u in v[-2::-1])
    return v


def xmatrix(g, format='vertical', include_basis_vectors=False):
    """
    Return the extra diagonals to add in each sector, as a matrix

    EXAMPLES::

        sage: xmatrix(2)
        [a 1]
        [a a]
        [1 a]
        sage: xmatrix(2, format='horizontal')
        [a a 1]
        [1 a a]
    """
    g, n, a, K, V = structure(g)
    b = basic_diagonals(g, all_the_way=include_basis_vectors)
    if format == 'vertical':
        return matrix(K, b[1:] + [u[::-1] for u in b[-2:0:-1]])
    elif format == 'horizontal':
        return matrix(K, b[1:] + [u[::-1] for u in b[-2:0:-1]]).transpose()
    raise ValueError("format should be 'horizontal' or 'vertical'")


def augmatrix(g, include_basis_vectors=False):
    """
    Return the augmentation matrix, ie diagonals to add in each sector

    Same as xmatrix, different implementation to compare speed.

    EXAMPLES::

        sage: augmatrix(2, include_basis_vectors=True)
        [1 a a 1 0]
        [0 1 a a 1]
    """
    g, n, a, K, V = structure(g)
    s = list(number_sequence(g, all_the_way=True))
    if include_basis_vectors:
        return matrix(K, (s[1:], s[:-1]))
    return matrix(K, (s[2:-1], s[1:-2]))


def within_bounds(v, xmax, ymax):
    return v[0] < xmax and v[1] < ymax


def long_diagonals(g, xmax=10, ymax=10, clip=True):
    """
    Return the long diagonals for 0 <= x < xmax, 0 <= y < ymax
    """
    g, n, a, K, V = structure(g)
    xmax, ymax = K(xmax), K(ymax)
    M = MatrixSpace(K, 2)
    def bounded(v):
        return within_bounds(v, xmax, ymax)
    x = identity_matrix(K, 2).rows()
    if xmax == ymax:
        x = basic_diagonals(g)
    m = xmatrix(g)
    l = m.nrows()
    k, nc = 1, len(x)
    while k < nc:
        if bounded(x[k-1]) and bounded(x[k]):
            xx = (m * M(x[k-1:k+1])).rows()
            if any(bounded(v) for v in xx):
                x[k:k] = xx
                nc += l
            else: k += 1
        else:
            k += 1
    if clip:
        x = list(filter(bounded, x))
    if xmax == ymax:
        x.extend(v[::-1] for v in x[-2::-1])
    return x


def low_diagonals(g=2, xmax=10, clip=True, start=None):
    """
    Return the long diagonals for 0 <= x < xmax and 0 <= y < x

    These are returned as a dictionary d whose keys are the
    x coordinates of these diagonals, and d[x] is the sorted
    list of y coordinates such that (x, y) is such a diagonal.

    If a starting list of (unclipped) diagonals is passed in,
    it is used as a starting point and extended in place.
    """
    g, n, a, K, V = structure(g)
    xmax = K(xmax)
    M = MatrixSpace(K, 2)
    def bounded(v):
        return v[0] <= xmax
    m = xmatrix(g)
    l = m.nrows()
    if start is None:
        x = basic_diagonals(g)
    else:
        x = start
    k, nc = 1, len(x)
    while k < nc:
        if bounded(x[k-1]) and bounded(x[k]):
            xx = (m * M(x[k-1:k+1])).rows()
            if any(bounded(v) for v in xx):
                x[k:k] = xx
                nc += l
            else: k += 1
        else:
            k += 1
    if clip:
        x = list(filter(bounded, x))
    d = defaultdict(list)
    for u, v in x:
        d[u].append(v)
    return d


def unclipped_very_low_diagonals(g=2, xmax=10, start=None, verbose=False):
    """
    Return the long diagonals for 0 <= x < xmax and 0 <= y < a/2 * x

    These are returned sorted by slope, unclipped (ie there are
    some with x > xmax, so that it can be used to start again.
    """
    g, n, a, K, V = structure(g)
    xmax = K(xmax)
    M = MatrixSpace(K, 2)
    def bounded(v):
        return v[0] <= xmax
    m = xmatrix(g)
    l = m.nrows()
    if start is None:
        x = basic_diagonals(g)
        k = len(x) - 1
        x[k:] = (m[:g] * M(x[k-1:k+1])).rows()
    else:
        x = start
    k, nc = 1, len(x)
    if verbose:
        indent = ' |' * verbose
        old_percent = 0
        old_time, old_iso = datestamp()
        start_time, start_iso = old_time, old_iso
        print(f'[{start_iso}] {indent}   0%')
    while k < nc:
        if verbose:
            pc = (200 * x[k][1] / a / x[k][0]).floor()
            if pc > old_percent:
                percent = pc
                new_time, new_iso = datestamp()
                this_pc = str(new_time - old_time).split('.')[0]
                all_pc = str(new_time - start_time).split('.')[0]
                timing = f'this % took {this_pc}, all % so far took {all_pc}'
                print(f'[{new_iso}]{indent} {percent:3d}% {timing}')
                old_percent = percent
                old_time, old_iso = new_time, new_iso
        if bounded(x[k-1]) and bounded(x[k]):
            xx = (m * M(x[k-1:k+1])).rows()
            if any(bounded(v) for v in xx):
                x[k:k] = xx
                nc += l
            else:
                k += 1
        else:
            k += 1
    return x


def low_diags_from_unclipped_very_low(vecs, xmax, only_very_low=False):
    r"""
    Only for genus two. May need more attention for general g.
    """
    V = vecs[0].parent()
    K = V.base_ring()
    a = K.gen()
    x = list(filter(lambda x: x[0] <= xmax, vecs))
    if not only_very_low:
        k = x.index(V((a, K.one())))
        x.extend((V((x, a*x-y)) for x, y in x[-2:k-1:-1]))
        assert x[-1][0] == x[-1][1]  # consistency check
    d = defaultdict(list)
    for u, v in x:
        d[u].append(v)
    return d


def length_ratios(g):
    """
    Return the length ratios of the diagonals to the long diagonals
    """
    ns = number_sequence(g)
    return [s - ns[i-1] for i, s in enumerate(ns) if i > 0]


def all_diagonals(g, xmax=10, ymax=10):
    kk = length_ratios(g)
    m = min(kk)
    d0 = long_diagonals(g, xmax=xmax/m, ymax=ymax/m, clip=True)
    d = []
    for k in kk:
        def bounded(v):
            return within_bounds(v, xmax/k, ymax/k)
        d.extend((k * v for v in filter(bounded, d0)))
    return d


def animate_long_diagonals(gmin=2, gmax=27, xymax=100, figsize=(10,10), verbose=True):
    """
    Animate long diagonals pictures for a range of ``g``.

    EXAMPLE::

        sage: a = animate_long_diagonals(gmin=2, gmax=100, xymax=300, figsize=(15, 15))  # long time
        sage: a  # long time
    """
    if verbose:
        import time
    diags = []
    if verbose:
        t = [time.time()]
    for g in sxrange(gmin, gmax+1):
        if verbose:
            print('[{datestamp()}] g = {}'.format(g))
        p = nf_point2d(long_diagonals(g, xymax, xymax),
            aspect_ratio=1, figsize=(10, 10))
        p.save('long_diagonals_g_{}_xy_{}.svg'.format(g,xymax))
        p.save('long_diagonals_g_{}_xy_{}.png'.format(g,xymax))
        diags.append(p)
        if verbose:
            t.append(time.time())
            print('[{datestamp()}]         took', int(t[-1] - t[-2]), 's')
    if verbose:
        print('g = {} to {} took {} s'.format(gmin, gmax, int(t[-1]-t[0])))
    a = animate(diags)
    a.save('long_diagonals_g_{}_to_{}_xy_{}.gif'.format(gmin, gmax, xymax))
    return a


def boxes(d, box_size=1):
    """
    Put points of the subset d of RR^2 into b by b square boxes

    EXAMPLES::

        sage: xmax = ymax = 2^9  # long time
        sage: b = 4  # long time
        sage: d2 = long_diagonals(2, xmax, ymax, clip=True)  # long time
        sage: len(d2)  # long time
        xxx
        sage: B = boxes(d2, box_size=b, return_box=True)  # long time
        sage: G = Graphics()  # long time
        sage: for s, l in box.iteritems():  # long time
        ....:     G += text('{0:02d}'.format(len(l)), s)
        sage: G.show(aspect_ratio=1, aspect_ratio=1, figsize=32)  # long time
        Launched png viewer for Graphics object consisting of ... graphics primitives
        sage: b = 4
        sage: from collections import defaultdict
        sage: stats = defaultdict(list)
        sage: for x, y in cartesian_product([range(0, xmax, b), range(0, ymax, b)]):  # long time
        ....:     stats[len(B[(x, y)])].append((x,y))
        sage: for n, l in stats.iteritems():  # long time
        ....:     print(n, len(l))
        0 2
        1 6
        2 40
        3 176
        4 403
        5 902
        6 1616
        7 2192
        8 2575
        9 2497
        10 2107
        11 1536
        12 1076
        13 650
        14 325
        15 148
        16 63
        17 52
        18 14
        19 2
        20 2

        sage: xmax = ymax = 2^10
        sage: b = 4
        sage: d2 = long_diagonals(2, xmax, ymax, clip=True)  # long time
        sage: len(d2)  # long time
        572605
        sage: B = boxes(d2, box_size=b, return_box=True)  # long time
        sage: G = Graphics()  # long time
        sage: for s, l in box.iteritems():  # long time
        ....:     G += text('{0:02d}'.format(len(l)), s)
        sage: G.show(aspect_ratio=1, aspect_ratio=1, figsize=64)  # long time
        Launched png viewer for Graphics object consisting of 65532 graphics primitives
        sage: 2^16
        65536
        sage: b = 4
        sage: from collections import defaultdict
        sage: stats = defaultdict(list)
        sage: for x, y in product(range(0, xmax, b), range(0, ymax, b)):  # long time
        ....:     stats[len(B[(x, y)])].append((x,y))
        sage: for n, l in stats.iteritems():  # long time
        ....:     print(n, len(l))
        ....:
        0 4
        1 48
        2 206
        3 772
        4 2034
        5 3812
        6 6413
        7 8486
        8 9836
        9 9623
        10 8368
        11 6248
        12 4158
        13 2684
        14 1489
        15 714
        16 351
        17 184
        18 72
        19 20
        20 10
        21 2
        22 2
        sage: from itertools import chain
        sage: for k, l in stats.iteritems():  # long time
        ....:     m = len(l)
        ....:     G = Graphics()
        ....:     x, y = vector(ZZ, l[randint(0, m-1)])
        ....:     xx = FiniteEnumeratedSet((x-b, x, x+b))
        ....:     yy = FiniteEnumeratedSet((y-b, y, y+b))
        ....:     xy = cartesian_product((xx, yy))
        ....:     G += nf_point2d(chain(B[s] for s in xy))
        ....:     G += line2d([(x, y), (x+4, y), (x+4, y+4), (x, y+4), (x, y)])
        ....:     G.show()
    """
    b = box_size
    from collections import defaultdict
    box = defaultdict(list)
    xmax = ymax = 2^4
    for x, y in d:
        box[(b*floor(x/b), b*floor(y/b))].append((x,y))
    return box


def nf_point2d(points, bits=20, **options):
    """
    Plot points with coordinates in a number field

    INPUT:

    - ``points`` -- a list or iterable of
      pairs of number field elements representing points

    - ``bits`` (default: 20): how many bits for the fractional part

    Usual options for :meth:`point2d` can be used, eg `color='green'`
    or `aspect_ratio=1`.

    Note that applying the usual point2d to points with coordinates
    in a number field can give extremely weird results, since the
    coordinates get converted to Python floats in a brutal way.
    When the coefficients for the representation of a number field
    element as a polynomial in the number field generator are huge,
    the error in converting it to a float can get way bigger than
    the actual number field element.
    """
    q = 2^bits
    qq = RDF(q)
    p = lambda x, y: (RDF((q*x).round())/qq, RDF((q*y).round())/qq)
    return point2d((p(x,y) for x, y in points), **options)


def RDF_from_nf(x, bits=20, **options):
    """
    Return an RDF elements from this number field element.

    INPUT:

    - ``l`` -- a list or iterable of number field elements

    - ``bits`` (default: 20): how many bits for the fractional part

    Usual options for :meth:`point2d` can be used, eg `color='green'`
    or `aspect_ratio=1`.
    """
    q = 2^bits
    qq = RDF(q)
    return RDF((q*x).round())/qq


def RDF_vector_from_nf_vector(u, bits=20, **options):
    """
    Return an RDF elements from this number field element.

    INPUT:

    - ``l`` -- a list or iterable of number field elements

    - ``bits`` (default: 20): how many bits for the fractional part

    Usual options for :meth:`point2d` can be used, eg `color='green'`
    or `aspect_ratio=1`.
    """
    q = 2^bits
    qq = RDF(q)
    return vector(RDF, ((q*x).round()/qq for x in u))


def RDF_list_from_nf_list(l, bits=20, **options):
    """
    Return a list of RDF elements from this list of number field elements.

    INPUT:

    - ``l`` -- a list or iterable of number field elements

    - ``bits`` (default: 20): how many bits for the fractional part
    """
    q = 2^bits
    qq = RDF(q)
    return [RDF_from_nf(x) for x in l]


def RDF_segment_list_from_nf_segment_list(s, bits=20, **options):
    """
    Return a list of RDF elements from this list of number field elements.

    INPUT:

    - ``s`` -- a list or iterable of 2-dimensional vectors over a number field

    - ``bits`` (default: 20): how many bits for the fractional part
    """
    q = 2^bits
    qq = RDF(q)
    return [(RDF_vector_from_nf_vector(x), RDF_vector_from_nf_vector(y))
            for (x, y) in s]
