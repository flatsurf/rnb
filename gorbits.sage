# -*- coding: utf-8 -*-

r"""
Double n-gon cutting sequences and regular n-gon periodic billiard paths

Enumerate periodic orbits for billiard in regular (2g+1)-gon,
coded as sequences of sides they hit.

Use transition diagrams, substitution, ...

We encode the substitutions that apply to the (2g) different
transition diagrams.

Each of them is stored as a dictionary, with the images of the
vertices 1 .. n, and the images of the relevant arrows ij.

.. TODO:

    - make classes for cuts words, path words, substitution words,
      that can eat strings
"""
from collections import namedtuple
from itertools import chain
from sage.geometry.polyhedron.constructor import Polyhedron

load('hecke_orbits.sage')

SegmentSequence = namedtuple("SegmentSequence", 'segments rotation_order cutting_word')

def update_globals(g):
    r"""
    Update global variables.
    """
    g, n, a, p, K, O, V, M = gnapKOVM(g)
    r = side_reflections(g)
    s = sector_substitutions(g)
    xm = augmatrix(g, include_basis_vectors=True)
    m = dict((k, xm[:, k:k+2]) for k in range(2*g))
    one_to_n = ''.join(str(k) for k in range(1, 2*g+2))
    _g = g
    _n = n
    _a = a
    _K = K
    _V = V
    _M = M
    _xm = xm
    _m = m
    _one_to_n = one_to_n
    _r = r
    _s = s
    return


@cached_function
def substitution_words(g):
    r"""
    Return the set of words over ``[0 .. 2*g - 1]``.

    EXAMPLES::

        sage: g = 2
        sage: substitution_words(g)
        Finite words on [0 .. 3]
    """
    W = FiniteWords(range(2*g))
    W.rename(f'Finite words on [0 .. {2*g-1}]')
    return W


@cached_function
def cutting_words(g):
    r"""
    Return the set of words over [1 .. 2*g].

    These are used to describe cutting sequences.

    EXAMPLES::

        sage: g = 2
        sage: C = cutting_words(g)
        sage: C
        Finite cutting sequences on [1 .. 5]
    """
    W = FiniteWords(2*g + 1)
    W.rename(f'Finite cutting sequences on [1 .. {2*g+1}]')
    return W


@cached_function
def path_words(g):
    r"""
    Return the set of path words over [1 .. 2*g].

    Here, we call path word a word consisting of
    node labels and arrow labels.

    Specifically, to a cutting word like (1, 2),
    we associate the path word (1, (1, 2), 2, (2, 1)).

    EXAMPLES::

        sage: g = 1
        sage: P = path_words(g)
        sage: P
        Path words on [1 .. 3]
        sage: P.alphabet()
        {1, 2, 3, word: 12, word: 21, word: 23, word: 32}
        sage: g = 2
        sage: P = path_words(g)
        sage: P
        Path words on [1 .. 5]
        sage: print(*(f'{a}' for a in P.alphabet()))
        1 2 3 4 5 12 21 23 32 34 43 45 54
    """
    W = cutting_words(g)
    one_to_n = W(list(W.alphabet()))
    forwards = list(one_to_n.factor_iterator(2))
    backwards = [w[::-1] for w in forwards]
    P = FiniteWords(list(W.alphabet()) + sorted(forwards + backwards))
    P.rename(f'Path words on [1 .. {2*g+1}]')
    return P


def fix_substitution_word(w, g=None, words=None):
    if words is None:
        if g is None:
            raise ValueError("Please specify g or words")
        return words(w)
    pass


def fix_cutting_dictionary(p, g=None):
    r"""
    Transform a dictionary (in place) to make it use words.

    This is to let other functions accept easier to input
    dictionaries, say using strings or tuples.

    EXAMPLE::

        sage: p = {'', '12'}
        sage: fix_path_dictionary(p, g=2)
        sage: p
        {word: : word: 12}
    """
    if g is None:
        raise ValueError("expected a positive integer for 'g'")
    if len(p) != 1:
        raise NotImplementedError("Only works for single-item dictionaries")
    (k, v), = p.items()
    if k not in ['', ()]:
        raise NotImplementedError(
            "Only works if dictionary key is empty string or empty tuple")
    S = substitution_words(g)
    C = cutting_words(g)
    a, b = v.split(',') if ',' in v else list(v)
    v = int(a), int(b)
    del p[k]
    p[S()] = C(v)


def zigzag_from_cyclic(g, s):
    r"""
    Return the zigzag label of the side with cyclic label ``s``.

    The double `(2g+1)`-gon has cyclic side labels ``0 .. (2g)``
    where oriented side ``i`` (cycling in counterclockwise order)
    makes angle ``i*tau/n`` with the horizontal. In particular
    side ``0`` is horizontal.

    The zigzag labeling, due to Diana Davis, zigzags as in
    Figure 21 of [DPU].

    Ref: DPU: Davis, Pasquinelli, Ulcigrai. Cutting sequences
    on Bouw-Möller surfaces: an S-adic characterization.
    Annales Scientifiques de l'École Normale Supérieure 52:4
    (2019), 927--1023. :arXiv:`1509.03905`.

    EXAMPLES::

        sage: g = 2
        sage: [zigzag_from_cyclic(g, i) for i in range(2*g+1)]
        [5, 4, 2, 1, 3]
        sage: g = 3
        sage: [zigzag_from_cyclic(g, i) for i in range(2*g+1)]
        [7, 6, 4, 2, 1, 3, 5]
    """
    n = 2*g + 1
    s = int(s)
    if s == 0:
        return n
    if s <= g:
        return 2*(g + 1 - s)
    if s >= 2*g + 1:
        raise ValueError(f"Expected 0 <= s <= 2*g; got s = {s}, g = {g}")
    return 2*(s - g) - 1


def cyclic_from_zigzag(g, s):
    r"""
    Return the cyclic label of the side with zigzag label ``s``.

    The double `(2g+1)`-gon has cyclic side labels ``0 .. (2g)``
    where oriented side ``i`` (cycling in counterclockwise order)
    makes angle ``i*tau/n`` with the horizontal. In particular
    side ``0`` is horizontal.

    The zigzag labeling, due to Diana Davis, zigzags as in
    Figure 21 of [DPU].

    Ref: DPU: Davis, Pasquinelli, Ulcigrai. Cutting sequences
    on Bouw-Möller surfaces: an S-adic characterization.
    Annales Scientifiques de l'École Normale Supérieure 52:4
    (2019), 927--1023. :arXiv:`1509.03905`.

    EXAMPLES::

        sage: g = 2
        sage: [cyclic_from_zigzag(g, str(i)) for i in range(1, 2*g+2)]
        [3, 2, 4, 1, 0]
        sage: g = 3
        sage: [cyclic_from_zigzag(g, str(i)) for i in range(1, 2*g+2)]
        [4, 3, 5, 2, 6, 1, 0]
    """
    s = int(s)
    if s < 1 or s > 2*g+1:
        raise ValueError(f"Expected 1 <= s <= 2*g + 1; got s = {s}, g = {g}")
    if s == 2*g + 1:
        return 0
    if s % 2:
        return g + (s + 1)//2
    return g + 1 - s//2


def side_reflections(g):
    r"""
    Return the side reflection dictionary.

    This dictionary records, for each side, a dictionary
    recording how every other side gets renumbered when
    the billiard bounces off that side.

    EXAMPLES::

        sage: arrows = lambda items: (f"{a}→{b}" for a, b in items)
        sage: s = side_reflections(2)
        sage: for k, d in sorted(s.items()):
        ....:     print(f'{k}:', ', '.join(arrows(sorted(d.items()))))
        1: 1→1, 2→3, 3→2, 4→5, 5→4
        2: 1→4, 2→2, 3→5, 4→1, 5→3
        3: 1→5, 2→4, 3→3, 4→2, 5→1
        4: 1→3, 2→5, 3→1, 4→4, 5→2
        5: 1→2, 2→1, 3→4, 4→3, 5→5
        sage: s = side_reflections(3)
        sage: for k, d in sorted(s.items()):
        ....:     print(f'{k}:', ', '.join(arrows(sorted(d.items()))))
        1: 1→1, 2→3, 3→2, 4→5, 5→4, 6→7, 7→6
        2: 1→4, 2→2, 3→6, 4→1, 5→7, 6→3, 7→5
        3: 1→5, 2→7, 3→3, 4→6, 5→1, 6→4, 7→2
        4: 1→7, 2→6, 3→5, 4→4, 5→3, 6→2, 7→1
        5: 1→6, 2→4, 3→7, 4→2, 5→5, 6→1, 7→3
        6: 1→3, 2→5, 3→1, 4→7, 5→2, 6→6, 7→4
        7: 1→2, 2→1, 3→4, 4→3, 5→6, 6→5, 7→7
    """
    res = dict()
    n = 2*g + 1
    zfc = lambda s: zigzag_from_cyclic(g, s)
    cfz = lambda s: cyclic_from_zigzag(g, s)
    for k in range(1, 2*g+2):
        res_k = dict()
        res_k[k] = k
        ck = cfz(k)
        for i in range(1, g + 1):
            j = zfc((ck - i) % n)
            jj = zfc((ck + i) % n)
            res_k[j] = jj
            res_k[jj] = j
        res[k] = res_k
    return res


def sector_substitutions(g):
    r"""
    Return the substitutions from transition diagrams.

    These are returned as a dictionary of dictionaries.

    We guess that a way to obtain them is as on the photo
    of the chalkboard in Kassar 204 on 2019-12-16.

    EXAMPLES:

    Check that we recover the double polygon substitutions::

        sage: s = sector_substitutions(2)
        sage: mappings = lambda items: (f"{a}→{b}" for a, b in items)
        sage: for k, d in s.items():
        ....:     letters = sorted((a, b) for a, b in d.items() if a in ZZ)
        ....:     arrows = sorted((a, b) for a, b in d.items() if a not in ZZ)
        ....:     print(f"{k}: {', '.join(mappings(letters))},\n"
        ....:           f"   {', '.join(mappings(arrows))}")
        0: 1→2, 2→1, 3→4, 4→3, 5→5,
           12→, 21→, 23→23, 32→32, 34→, 43→, 45→4, 54→4
        1: 1→3, 2→5, 3→1, 4→4, 5→2,
           12→4, 21→4, 23→432, 32→234, 34→23, 43→32, 45→3, 54→3
        2: 1→4, 2→2, 3→5, 4→1, 5→3,
           12→3, 21→3, 23→34, 32→43, 34→432, 43→234, 45→2, 54→2
        3: 1→1, 2→3, 3→2, 4→5, 5→4,
           12→2, 21→2, 23→, 32→, 34→34, 43→43, 45→, 54→
        sage: s = sector_substitutions(3)
        sage: mappings = lambda items: (f"{a}→{b}" for a, b in items)
        sage: for k, d in s.items():
        ....:     letters = sorted((a, b) for a, b in d.items() if a in ZZ)
        ....:     arrows = sorted((a, b) for a, b in d.items() if a not in ZZ)
        ....:     print(f"{k}: {', '.join(mappings(letters))},\n"
        ....:           f"   {', '.join(mappings(arrows))}")
        0: 1→2, 2→1, 3→4, 4→3, 5→6, 6→5, 7→7,
           12→, 21→, 23→23, 32→32, 34→, 43→, 45→45,
           54→54, 56→, 65→, 67→6, 76→6
        1: 1→5, 2→7, 3→3, 4→6, 5→1, 6→4, 7→2,
           12→6, 21→6, 23→654, 32→456, 34→45, 43→54, 45→5432,
           54→2345, 56→23, 65→32, 67→3, 76→3
        2: 1→4, 2→2, 3→6, 4→1, 5→7, 6→3, 7→5,
           12→3, 21→3, 23→345, 32→543, 34→5432, 43→2345, 45→23456,
           54→65432, 56→654, 65→456, 67→4, 76→4
        3: 1→3, 2→5, 3→1, 4→7, 5→2, 6→6, 7→4,
           12→4, 21→4, 23→432, 32→234, 34→23456, 43→65432, 45→6543,
           54→3456, 56→345, 65→543, 67→5, 76→5
        4: 1→6, 2→4, 3→7, 4→2, 5→5, 6→1, 7→3,
           12→5, 21→5, 23→56, 32→65, 34→6543, 43→3456, 45→34,
           54→43, 56→432, 65→234, 67→2, 76→2
        5: 1→1, 2→3, 3→2, 4→5, 5→4, 6→7, 7→6,
           12→2, 21→2, 23→, 32→, 34→34, 43→43, 45→,
           54→, 56→56, 65→65, 67→, 76→
    """
    zfc = lambda s: zigzag_from_cyclic(g, s)
    cfz = lambda s: cyclic_from_zigzag(g, s)
    n = 2*g + 1
    res = dict()
    C = cutting_words(g)
    for k in range(2*g):
        res_k = dict()
        cfix_k = (k * g) % n
        # node images
        for i in range(-g, g + 1):
            res_k[zfc(i % n)] = zfc((cfix_k - i) % n)
        # edge images
        for i in range(1, n):
            j = res_k[i]
            jj = res_k[i+1]
            t = C()
            if jj - j > 1:
                t = C(range(j + 1, jj))
            if j - jj > 1:
                t = C(range(j - 1, jj, -1))
            tt = t[::-1]
            res_k[C([i, i+1])] = t
            res_k[C([i+1, i])] = tt
        res[k] = res_k
    return res


def side_reflections(g):
    r"""
    Return side reflections for the (2*g+1)-gon.

    Sides are numbered in zigzag numbering.

    EXAMPLES::

        sage: r = side_reflections(2)
        sage: arrows = lambda items: (f"{a}→{b}" for a, b in items)
        sage: for k, d in sorted(r.items()):
        ....:     print(f'{k}:', ', '.join(arrows(sorted(d.items()))))
        1: 1→1, 2→3, 3→2, 4→5, 5→4
        2: 1→4, 2→2, 3→5, 4→1, 5→3
        3: 1→5, 2→4, 3→3, 4→2, 5→1
        4: 1→3, 2→5, 3→1, 4→4, 5→2
        5: 1→2, 2→1, 3→4, 4→3, 5→5
        sage: r = side_reflections(3)
        sage: for k, d in sorted(r.items()):
        ....:     print(f'{k}:', ', '.join(arrows(sorted(d.items()))))
        1: 1→1, 2→3, 3→2, 4→5, 5→4, 6→7, 7→6
        2: 1→4, 2→2, 3→6, 4→1, 5→7, 6→3, 7→5
        3: 1→5, 2→7, 3→3, 4→6, 5→1, 6→4, 7→2
        4: 1→7, 2→6, 3→5, 4→4, 5→3, 6→2, 7→1
        5: 1→6, 2→4, 3→7, 4→2, 5→5, 6→1, 7→3
        6: 1→3, 2→5, 3→1, 4→7, 5→2, 6→6, 7→4
        7: 1→2, 2→1, 3→4, 4→3, 5→6, 6→5, 7→7
    """
    n = 2*g + 1
    res = dict()
    zfc = lambda s: zigzag_from_cyclic(g, s)
    cfz = lambda s: cyclic_from_zigzag(g, s)
    Zmodn = Zmod(n)
    for k in Zmodn:
        res_k = dict()
        for j in Zmodn:
            res_k[zfc(j)] = zfc(Zmodn(2*k-j))
        res[zfc(k)] = res_k
    return res


def apply_subs(i, zw, g='genus', s=None):
    """
    Return the result of applying substitution number i to this cutting word.

    INPUT:

    - ``i`` -- the index of the substitution to apply

    - ``zw`` -- the zigzag word to which to apply the substitution

    - ``g`` -- the genus

    - ``s`` -- a list or tuple of substitutions

    EXAMPLES:

    Apply all substitutions to word (1, 2) in varying genus::

        sage: # load('gorbits.sage')
        sage: for g in range(6):
        ....:     print(*(f'{apply_subs(k, (1, 2), g=g, s=None)}' for k in range(2*g)))
        21 1232
        21 3454 4323 1232
        21 5676 4323 3454 6545 1232
        21 7898 4323 5676 6545 3454 8767 1232
        21 9,10,11,10 4323 7898 6545 5676 8767 3454 10,9,8,9 1232

    Apply substitution ``0`` to word (1, 2) in varying genus::

        sage: print(*(f'{apply_subs(0, (1, 2), g=k, s=None)}' for k in range(1, 7)))
        21 21 21 21 21 21

    Apply substitution ``1`` to word (1, 2) in varying genus::

        sage: print(*(f'{apply_subs(1, (1, 2), g=k, s=None)}' for k in range(1, 7)))
        1232 3454 5676 7898 9,10,11,10 11,12,13,12
    """
    if s is None:
        s = sector_substitutions(g)
    subs = s[i]
    C = cutting_words(g)
    n = len(zw)
    zwz = C(chain(zw, zw[:1]))
    arrows = [zwz[k:k+2] for k in range(n)]
    return prod((C([subs[u]]) * C(subs[arrows[k]])
                 for k, u in enumerate(zw)), C())


def augment(p, w, g='genus', s=None):
    r"""
    Augment the collection of descendents of a periodic cutting sequence.

    INPUT:

    - ``p`` -- dictionary of periodic cutting sequence words

    - ``w`` -- word in the substitutions

    - ``s`` -- dictionary (default: ``_s``) -- substitution tuple

    EXAMPLES:

    Apply substitution words [0] and [0, 1] to zigzag word [1, 2]::

        sage: # load('gorbits.sage')
        sage: for g in range(2, 8):
        ....:     Z = cutting_words(g)
        ....:     S = substitution_words(g)
        ....:     p = p12 = {S(): Z([1, 2])}
        ....:     augment(p, S([0, 1]), g=g)
        ....:     augment(p, S([1, 0]), g=g)
        ....:     print(p[S([0])], p[S([1])], p[S([1, 0])], )
        21 3454 434543
        21 5676 656765
        21 7898 878987
        21 9,10,11,10 10,9,10,11,10,9
        21 11,12,13,12 12,11,12,13,12,11
        21 13,14,15,14 14,13,14,15,14,13

    Apply a long substitution word to zigzag words [1, 2] and [3, 4]::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: Z = cutting_words(g)
        sage: S = substitution_words(g)
        sage: w = S([1, 2, 3, 1, 3, 2, 2, 0, 0, 3, 0, 2, 2])
        sage: p = p12 = {S(): Z([1, 2])}
        sage: pp = p34 = {S(): Z([3, 4])}
        sage: augment(p, w, g=2, s=None)
        sage: augment(pp, w, g=2, s=None)
        sage: a, b, c = len(w), len(p), len(p[w])
        sage: d, e, f = len(w), len(pp), len(pp[w])
        sage: print(a, b, c, d, e, f)
        13 14 106628 13 14 172530
        sage: continued_fraction(f/c)
        [1; 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 184, 2]

    Apply all substitution words of length up to 4::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: Z = cutting_words(g)
        sage: S = substitution_words(g)
        sage: p = p12 = {S(): Z([1, 2])}
        sage: for w in S.iterate_by_length(4):
        ....:     augment(p, w, g=g, s=None)
        sage: len(p)
        341
        sage: sum(4**k for k in range(5))
        341
        sage: len(p[S([1, 3, 2, 1])])
        64
        sage: print(p[S([1, 3, 2, 1])])
        5432123432343212345434543212345434543212345434543212343234321234

    And now up to length 8::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: Z = cutting_words(g)
        sage: S = substitution_words(g)
        sage: p = p12 = {S(): Z([1, 2])}
        sage: pp = p34 = {S(): Z([3, 4])}
        sage: for w in S.iterate_by_length(8):  # long time
        ....:     augment(p, w, g=g, s=None)
        ....:     augment(pp, w, g=2, s=None)
        sage: len_12 = list(map(len, p12.values()))  # long time
        sage: len_34 = list(map(len, p34.values()))  # long time
        sage: print(len(p12), len(p34), max(len_12), max(len_34))  # long time
        87381 87381 9160 14822
        sage: len_12 = sorted(len_12)  # long time
        sage: len_34 = sorted(len_34)  # long time
        sage: len_all = sorted(len_12 + len_34)  # long time
        sage: len(len_all)  # long time
        174762

    Some histograms::

        sage: histogram(len_12, bins=4580)  # long time
        Graphics object consisting of 1 graphics primitive
        sage: histogram(len_12, bins=4580/32)  # long time
        Graphics object consisting of 1 graphics primitive
        sage: histogram(len_34, bins=7411)  # long time
        Graphics object consisting of 1 graphics primitive
        sage: histogram(len_34, bins=7411/32)  # long time
        Graphics object consisting of 1 graphics primitive
        sage: histogram(len_all, bins=7411)  # long time
        Graphics object consisting of 1 graphics primitive
        sage: histogram(len_all, bins=7411/32)  # long time
        Graphics object consisting of 1 graphics primitive
    """
    if len(p) == 1 and list(p.keys())[0] in ['', ()]:
        fix_cutting_dictionary(p, g=g)
    if s is None:
        s = sector_substitutions(g)
    if w not in p:
        if w[:-1] not in p:
            augment(p, w[:-1], g=g, s=s)
        if max(p[w[:-1]].parent().alphabet()) != 2*g+1:
            raise ValueError("Augmenting dictionary of the wrong genus")
        p[w] = apply_subs(w[-1], p[w[:-1]], g=g, s=s)


def fold_once(u, v, g='genus', r=None):
    r"""
    Fold a cutting sequence once.

    INPUT:

    - ``u`` -- part of the path already done

    - ``v`` -- period of cutting sequence

    - ``g`` -- genus

    - ``r`` -- dictionary of side reflections

    The trajectory hits the first side in ``v``, so
    - the first letter of ``v`` gets appended to ``u``,
    - the new ``v`` has its first letter shifted to the end,
      and the reflection corresponding to hitting that side
      is applied to all letters of ``v``.

    .. NOTE:

        This function is nowhere used. Can probably be erased.

    EXAMPLE::

        sage: g = 2
        sage: C = cutting_words(g)
        sage: fold_once(C(), C([1, 2]), g=g)
        (word: 1, word: 31)
    """
    if r is None:
        r = side_reflections(g)
    C = cutting_words(g)
    a = v[0]
    ra = WordMorphism(r[a], domain=C, codomain=C)
    return u * C([a]), ra(v[1:]) * v[:1]


def fold(u, g='genus', r=None):
    r"""
    Fold double pentagon cutting sequence into pentagon billiard orbit.

    This is done by keeping a current letter and
    transformation.

    Could maybe be made faster by stopping at half the length
    when the path has the five-fold symmetry. (Maybe not worth it.)

    INPUT:

    - ``u`` -- a cutting sequence word

    - ``g`` -- genus

    - ``r`` -- reflection dictionary

    OUTPUT:

    Folded up word.

    EXAMPLES::

        sage: print(fold((1, 2), g=2, r=None))
        1245312453
        sage: print(fold((3, 4), g=2, r=None))
        sage: print(fold((1, 2), g=5, r=None))
        1,2,4,6,8,10,11,9,7,5,3,1,2,4,6,8,10,11,9,7,5,3
        sage: print(fold((3, 4), g=2, r=None))
        3415234152
        sage: print(fold((1, 2), g=4, r=None))
        124689753124689753
        sage: print(fold((3, 4), g=4, r=None))
        349349
        sage: print(fold((5, 6), g=4, r=None))
        563819274563819274
        sage: print(fold((7, 8), g=4, r=None))
        784159623784159623
        sage: print(fold((7, 8), g=7, r=None))  # 3-fold symmetry
        7,8,9,6,11,4,13,2,15,1,14,3,12,5,10,7,8,9,6,11,4,13,2,15,1,14,3,12,5,10
        sage: print(fold((9, 10), g=7, r=None))  # 5-fold symmetry
        9,10,3,15,4,9,10,3,15,4
    """
    if r is None:
        r = side_reflections(g)
    C = cutting_words(g)
    u = C(u)
    a, b = u[:1], u[1:]
    x = b[:1]
    rr = C.identity_morphism()
    r0 = WordMorphism(r[a[0]], domain=C, codomain=C)
    while not (rr == r0 and b == u):
        a *= x
        rr = WordMorphism(r[x[0]], domain=C, codomain=C) * rr
        b = b[1:]
        if not b:
            b = u[:]
        x = rr(b[:1])
    return a


def vector_of_word(w, g='genus', v=None, V=None, m=None):
    """
    Return vector associated to substitution word.

    Start with vector (1, 0). Apply sector matrices
    following the given word. Get direction vector
    in staircase picture.

    In other words this is the vector for that word
    that is in the orbit of (1, 0) under the
    Veech group of the staircase.

    INPUT:

    - ``w`` -- a substitution word

    - ``g`` -- the genus

    - ``v`` -- a starting vector

    - ``V`` -- vector space

    - ``m`` -- dictionary of sector matrices

    If ``v`` is ``None``, start with vector ``(1, 0)``.

    EXAMPLE::

        sage: w = (0, 0)
        sage: vector_of_word(w, g=2)
        (1, 0)
        sage: [vector_of_word(x, g=2) for x in ([0], [1], [2], [3])]
        [(1, 0), (a, 1), (a, a), (1, a)]

    .. TODO:

        - Add example with starting vector other than (1, 0)
    """
    if V is None:
        g, n, a, p, K, O, V, M = gnapKOVM(g)
    if m is None:
        xm = augmatrix(g, include_basis_vectors=True)
        m = dict((k, xm[:, k:k+2]) for k in range(2*g))
    if not v:
        v = (1, 0)
    v = V(v)
    for x in w:
        v = m[x] * v
    return v


def which_sector(v, g='genus', xm=None):
    r"""
    Return the sector number for this vector.

    INPUT:

    - ``v`` -- vector

    - ``g`` -- genus

    - ``xm`` -- augmentation matrix

    EXAMPLES::

        sage: g, n, a, K, V = structure(2)
        sage: vectors = [V([2, 0]), V([2, 1]), V([3, 2]), V([2, 2]), V([2, 3]), V([1, 3])]
        sage: [which_sector(V(v), g=2) for v in vectors]
        [0, 0, 1, 2, 2, 3]
    """
    if xm is None:
        xm = augmatrix(g, include_basis_vectors=True)
        xm[:,1:1] = xmatrix(g, format='horizontal')
    x, y = v
    if y == 0:
        return 0
    if x == 0:
        raise ValueError("Vertical vector is out of the cone.")
    sector = 0
    if x > 0:
        slope = y/x
        while sector < 2*g - 1 and slope >= xm[1, sector + 1] / xm[0, sector + 1]:
            sector += 1
        return sector


def stretch_sector(v, g='genus', m=None):
    r"""
    Apply the inverse matrix of the sector matrix
    for this vector to this vector.

    EXAMPLES::

        sage: g = 2
        sage: g, n, a, p, K, O, V, M = gnapKOVM(g)
        sage: stretch_sector(v=V([3,2]), g=g)
        (a, 2*a - 3)
        sage: g = 5
        sage: g, n, a, K, V = structure(g)
        sage: stretch_sector(V([2,2]), g=g)
        (2*a^4 - 2*a^3 - 6*a^2 + 4*a + 2, 0)
    """
    g, n, a, p, K, O, V, M = gnapKOVM(g)
    if m is None:
        xm = augmatrix(g, include_basis_vectors=True)
        m = dict((k, xm[:, k:k+2]) for k in range(2*g))
    v = V(v)
    return (m[which_sector(v, g=g)])^-1 * v


def word_of_vector(v, g='genus'):
    r"""
    Return the substitution word that maps (1, 0) to this vector.

    EXAMPLE::

        sage: g = 2
        sage: g, n, a, p, K, O, V, M = gnapKOVM(g)
        sage: vectors = [(1, 0), (a, 1), (a, a), (1, a)]
        sage: print(*(f"'{word_of_vector(v, g=g)}'" for v in vectors))
        '' '1' '2' '3'
    """
    vec = v
    S = substitution_words(g)
    itinerary = []
    while vec[1] > 0:
        itinerary.append(which_sector(vec, g=g))
        vec = stretch_sector(vec, g=g)
    return S(itinerary[::-1])

def sectors_of_vector(v, g='genus'):
    r"""
    Iterator of substitution words from this vector.

    EXAMPLE::

        sage: g = 2
        sage: g, n, a, p, K, O, V, M = gnapKOVM(g)
        sage: uu = [(1, 0), (a, 1), (a, a), (1, a)]
        sage: v, w, x, y = [sectors_of_vectors(u) for u in uu]
        sage: [next(w), next(x), next(y)]
        [1, 2, 3]
        sage: next(v)
        Traceback (most recent call last):
        ...
        StopIteration...
    """
    vec = v
    S = substitution_words(g)
    while vec[1] > 0:
        s = which_sector(vec, g=g)
        vec = stretch_sector(vec, g=g)
        yield s


def midpoint_reverse_zigzag(g):
    r"""
    Return midpoint coordinates in reverse zigzag order.

    The origin is at the midpoint of the longest side
    (the diagonal of the square) in the staircase polygon.

    EXAMPLES::

        sage: midpoint_reverse_zigzag(2)
        [(a + 1, -a), (a + 1, 0), (0, 0), (0, a + 1), (-a, a + 1)]
    """
    aa = number_sequence(g, all_the_way=True)
    zero = aa[0]
    V = VectorSpace(zero.parent(), 2)
    mx = [-sum(aa[1:k:2], aa[0]) - sum(aa[3:k+2:2], aa[0]) for k in range(2*g + 1)]
    my = [sum(aa[0:k:2], aa[0]) + sum(aa[2:k+2:2], aa[0]) for k in range(2*g + 1)]
    mm = [V(v) for v in zip(mx, my)]
    o = mm[g]
    mm = [v - o for v in mm]
    return mm


def midpoint_translate_side(g, labels='zigzag'):
    r"""
    Return midpoints, translation vectors, and side equations in staircase view.

    The output consists in three dictionaries midpoint, translate, side
    whose keys are the side labels in zigzag or cyclic numbering, and
    values are: midpoints of enter sides as vectors over the underlying
    number field, translation vectors from exit to enter side as vectors
    over the underlying number field, and equations for the lines supporting
    the exit sides, given as triples (a, b, c) representing the line equation
    `a*x + b*y = c` in the plane.

    The scale is such that the shortest side of the staircase polygon has length 2.

    EXAMPLES::

        sage: m, t, s = midpoint_translate_side(2, labels='cyclic')
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(m.items())))
        0: (a + 1, -a)
        1: (-a - 1, 0)
        2: (0, -a - 1)
        3: (-a, a + 1)
        4: (0, 0)
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(t.items())))
        0: (2*a + 2, -2*a)
        1: (-2*a - 2, 0)
        2: (0, -2*a - 2)
        3: (-2*a, 2*a + 2)
        4: (0, 0)
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(s.items())))
        0: (0, 1, a)
        1: (a, 1, 2*a + 1)
        2: (a, a, 0)
        3: (1, a, 2*a + 1)
        4: (1, 0, a)
        sage: m, t, s = midpoint_translate_side(2, labels='zigzag')
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(m.items())))
        1: (-a, a + 1)
        2: (0, -a - 1)
        3: (0, 0)
        4: (-a - 1, 0)
        5: (a + 1, -a)
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(t.items())))
        1: (-2*a, 2*a + 2)
        2: (0, -2*a - 2)
        3: (0, 0)
        4: (-2*a - 2, 0)
        5: (2*a + 2, -2*a)
        sage: print('\n'.join(f'{k}: {v}' for k, v in sorted(s.items())))
        1: (1, 0, a)
        2: (1, a, 2*a + 1)
        3: (a, a, 0)
        4: (a, 1, 2*a + 1)
        5: (0, 1, a)
    """
    # The side lengths are a_0 = 0, a_1 = 1, a_2 = a, ..., a_n = 0
    aa = number_sequence(g, all_the_way=True)
    zero = aa[0]
    V = VectorSpace(zero.parent(), 2)
    mm = midpoint_reverse_zigzag(g)
    mm = [-v if k % 2 else v for k, v in enumerate(mm)]
    # side equations
    n = 2*g + 1
    side_normal = [V(aa[k:k+2]) for k in range(n)]
    for k in range(1, n, 2):
        side_normal[k] = side_normal[k][::-1]
    side_equation = [(v[0], v[1], -v*mm[k])
                     for k, v in enumerate(side_normal)]
    if labels == 'zigzag':
        midpoint = dict((n - k, m) for k, m in enumerate(mm))
        translate = dict((k, 2*v) for k, v in midpoint.items())
        side = dict((n - k, m) for k, m in enumerate(side_equation))
        return midpoint, translate, side
    if labels == 'cyclic':
        midpoint = dict(enumerate(mm[:1] + mm[1::2] + mm[:0:-2]))
        translate = dict((k, 2*v) for k, v in midpoint.items())
        side = dict(enumerate(side_equation))
        return midpoint, translate, side
    raise ValueError(f"labels should be 'cyclic' or 'zigzag', got: '{labels}'")
    return


def staircase_ngon_vertices(g, labels=None):
    r"""
    Return coordinates of staircase ngon vertices.

    EXAMPLE::

        sage: staircase_ngon_vertices(2)
        [(a, -a), (a + 2, -a), (a, a), (-a, a + 2), (-a, a)]
        sage: staircase_ngon_vertices(2, labels='cyclic')
        {0: (a, -a), 1: (a + 2, -a), 2: (a, a), 3: (-a, a + 2), 4: (-a, a)}
        sage: vv = staircase_ngon_vertices(2, labels='zigzag').items()
        sage: print(', '.join(f'{k}: {v}' for k, v in sorted(vv)))
        1: (-a, a + 2), 2: (a, a), 3: (-a, a), 4: (a + 2, -a), 5: (a, -a)
    """
    mm = midpoint_reverse_zigzag(g)
    V = mm[0].parent()
    ngon = [mm[0] - V((1, 0))]
    for m in mm[:1] + mm[1::2] + mm[:2:-2]:
        ngon.append(2*m - ngon[-1])
    if labels == 'cyclic':
        return dict(enumerate(ngon))
    if labels == 'zigzag':
        zfc = lambda k: zigzag_from_cyclic(g, k)
        return dict((zfc(k), v) for k, v in enumerate(ngon))
    return ngon


def hit_from_point_vector_to_line(p, d, e, V=None):
    r"""
    Return the hit we get starting from a point
    along a direction and getting to a line.

    INPUT:

    - ``p`` -- starting point

    - ``v`` -- vector for the direction of travel

    - ``e`` -- triple (a, b, c) representing the line
      of equation ``a * x + b * y == c``.

    METHOD:

    From a point ``p = (s, t)`` in the plane, travel in direction
    ``d = (u, v)``. Since a point ``(x, y)`` is on this line iff
    the displacement vector ``(x - s, y - t)`` is collinear to the
    direction vector ``d = (u, v)``, which happens iff the
    determinant of these two vectors is zero, the line equation
    for the trajectory is ``A * x + B * y == C``
    with ``A == v``, ``B == -u``, ``C == s*v - u*t``.

    The lines of equations ``A * x + B * y == C`` and
    ``a * x + b * y == c`` intersect at the point whose
    coordinates ``(x, y)`` are given by Cramer's formulas:
    ``x == (C*b - c*B) / (A*b - a*B)``,
    ``y == (A*c - c*A) / (A*b - a*B)``.

    EXAMPLES:

    We test where the line from point ``(0, 2)`` in
    direction ``[1, -1]`` intersects the line ``x - y == 0``::

        sage: p = vector((0,2))
        sage: d = vector((1,-1))
        sage: e = vector((1,-1,0))
        sage: hit_from_point_vector_to_line(p, d, e)
        (1, 1)

    We test where the line from point ``(0, -2)`` in
    direction ``[1, 2]`` intersects the line ``y == 2``::

       sage: p = vector((0,-2))
       sage: d = vector((1,2))
       sage: e = vector((0,1,2))
       sage: hit_from_point_vector_to_line(p, d, e)
       (2, 2)
    """
    s, t = p
    u, v = d
    a, b, c = e
    A, B, C = v, -u, v*s - u*t
    delta = (A * b - a * B)
    hit = (C*b - c*B) / delta, (A*c - a*C) / delta
    if V:
        return V(hit)
    return hit


def dp_orbit_segments(p, w, g='genus', v=None, verbose=False):
    """
    Orbit segments in staircase ngon for substitution word ``w``.

    INPUT:

    - ``p`` -- a dictionary mapping "contraction words" to "cutting sequence words"

    - ``w`` -- a "contraction word"

    - ``v`` (default: `None`) -- ...

    OUTPUT:

    - A list of orbit segments (as pairs of points)

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: G = Graphics()
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: seg = dp_orbit_segments(p, S([2, 1, 0]), g=g)
        sage: pict = sum((line2d(s) for s in seg), Graphics())
        sage: pict.show(aspect_ratio=1, figsize=10)
    """
    midpoint, translate, side = midpoint_translate_side(g, labels='zigzag')
    V = midpoint[g].parent()
    exit = lambda p, d, i: hit_from_point_vector_to_line(p, d, side[i], V=V)
    S = substitution_words(g)
    w = S(list(map(Integer, w.split(',') if ',' in w else list(w))))
    if w not in p:
        # cutting sequence for word w not computed yet; compute it
        augment(p, w, g=g)
    if verbose:
        print(f'w = {w}')
    s = p[w]  # the cutting sequence
    if verbose:
        print("%s: %s" %(w, p[w]))
    v = vector_of_word(w, g=g, v=v)
    sp = midpoint[s[0]]
    if verbose:
        print(f"starting point: {sp}")
    seg = []
    for i in s[1:]:
        ep = exit(sp, v, i)
        seg.append((sp, ep))
        sp = ep + translate[i]
    # last step: close up path to very first start point
    ep = midpoint[s[0]] - translate[s[0]]
    seg.append((sp, ep))
    return seg


def draw_dp_orbit(p, w, g='genus', v=None,
                  verbose=False,
                  line_thickness=1,
                  adapt_line_thickness=True,
                  ngon_thickness=1,
                  line_color='black',
                  ngon_color='blue'):
    """
    Return picture of double regular pentagon with orbit plotted

    The picture is returned as a Sage graphics object with
    aspect ratio 1 and no axes.

    One can set the colors of the n-gon outline and the line segments.

    EXAMPLES:

    Use a dictionary with (substitution word, cutting word)
    as key-value pairs::

        sage: # load('gorbits.sage')
        sage: g = 5
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: draw_dp_orbit(p, S([1, 2]), g=g, v=None, verbose=False, line_thickness=1)
        Graphics object consisting of 14 graphics primitives

    Or simply use strings::

        sage: # load('gorbits.sage')
        sage: g = 5
        sage: p = p12 = {'': '1, 2')}
        sage: draw_dp_orbit(p, '1, 0, 2', g=g, v=None, verbose=False, line_thickness=1)
        Graphics object consisting of 14 graphics primitives
    """
    seg = dp_orbit_segments(p, w, g=g, v=v, verbose=verbose)
    aa = number_sequence(g, all_the_way=True)
    V = VectorSpace(aa[0].parent(), 2)
    n = 2*g + 1
    mm = midpoint_reverse_zigzag(g)
    ngon = [mm[0] - V((1, 0))]
    for m in mm[:1] + mm[1::2] + mm[:0:-2]:
        ngon.append(2*m - ngon[-1])
    pi_over_n = RDF.pi()/n
    M = matrix(RDF, 2, [1, pi_over_n.cos(), 0, pi_over_n.sin()])
    ngon = [M * V(x) for x in ngon]
    rgon = [-v for v in ngon]
    GG = Graphics()
    options = {'aspect_ratio': 1, 'axes': False, 'color': ngon_color, 'thickness': ngon_thickness}
    GG += line2d(ngon, **options)
    GG += line2d(rgon, **options)
    if adapt_line_thickness:
        orbit_length = sum(norm(y - x) for x, y in seg)/norm(ngon[0])
        # thickness_from_length = lambda x: 0.05 + 0.95/(1 + (x/100)^2)
        thickness_from_length = lambda x: 2/(1 + x/100)
        line_thickness = thickness_from_length(orbit_length)
    GG = sum((line2d([M * x, M * y], thickness=line_thickness, color=line_color)
              for x, y in seg), GG)
    # GG += line2d(map(lambda x: M * x, [V(x) for x in mm]), color='red')
    return GG


def reflect(g):
    r"""
    Return the matrices for (affine) reflections in the
    sides of the staircase n-gon.

    This is returned as a dictionary whose keys are the
    zigzag labels for the n-gon sides.

    These affine reflections are the conjugates of the
    orthogonal reflections in the regular n-gon picture.

    EXAMPLE::

        sage: for k, m in sorted(reflect(2).items()):
        ....:     print(f'{k}:')
        ....:     print(f'{m}')
        1:
        [-1  0]
        [ a  1]
        2:
        [ a  1]
        [-a -a]
        3:
        [ 0 -1]
        [-1  0]
        4:
        [-a -a]
        [ 1  a]
        5:
        [ 1  a]
        [ 0 -1]
    """
    res = dict()
    vv = staircase_ngon_vertices(g, labels='cyclic')
    o = sum(vv.values())/len(vv)
    zfc = lambda k: zigzag_from_cyclic(g, k)
    n = 2*g + 1
    for k in range(n):
        a = vv[k] - o
        b = vv[(k + 1) % n] - o
        m = matrix([a, b])
        res[zfc(k)] = -m.transpose() * m[::-1].transpose().inverse()
    return res


def ngon_billiard_segments(p, w, g='genus', v=None,
                           print_cutting_word=False,
                           print_length=False,
                           distance_along=None,
                           return_rotation_order=False,
                           return_cutting_word=False):
    """
    Orbit segements in slanted n-gon billiard for substitution word ``w``.

    INPUT:

    - ``p`` -- a path dictionary with keys substitution words and values cutting words

    - ``w`` -- a substitution word

    - ``g`` -- genus (the ngon has `n = 2 g + 1`)

    - ``v`` -- initial vector (default: ``None``)

    - ``print_cutting_word`` -- whether to print the cutting word (default: ``False``)

    - ``print_length`` -- whether to print length (default: ``False``)

    - ``distance_along`` -- how far along the starting side to start,
      as a proportion of the side length, from 0 to 1 (default: ``None`` -- in that case
      we start from the midpoint of a side)

    - ``return_rotation_order`` -- whether to return the rotation order (default: ``False``)

    - ``return_billiard_word`` -- whether to return the billiard word (default: ``False``)

    OUTPUT: a named tuple with the following components:

    - ``segments`` -- the billiard orbit segments

    - ``rotation_order`` -- the rotation order, if requested
      (1 for no rotational symmetry, 3 for 3-fold rotational symmetry, etc.)

    - ``billiard_word`` -- the billiard word, if requested

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: seg = dp_orbit_segments(p, S([2, 1, 0]), g=g)
        sage: pict = sum((line2d(s) for s in seg), Graphics())
        sage: pict.show(aspect_ratio=1, figsize=10)
        sage: G = Graphics()
        sage: W = VectorSpace(RDF, 2)
        sage: rseg = [(W(a), W(b)) for (a, b) in seg]
        sage: for x in rseg:
        ....:     G += line2d(x)
        sage: ngon = staircase_ngon_vertices(g, labels='zigzag').values()
        sage: G += polygon(ngon, fill=False, color='red')
        sage: G.show(aspect_ratio=1, axes=False)
        sage: pi_over_n = RDF.pi()/n
        sage: M = matrix(RDF, 2, [1, pi_over_n.cos(), 0, pi_over_n.sin()])
        sage: GG = sum((line2d([M * x, M * y]) for x, y in seg),
        ....:          polygon((M * x for x in ngon), fill=False, color='red'))
        sage: GG.show(aspect_ratio=1, axes=False)
    """
    g, n, _, _, _, _, V, _ = gnapKOVM(g)
    zfc = lambda k: zigzag_from_cyclic(g, k)
    cfz = lambda k: cyclic_from_zigzag(g, k)
    midpoint, translate, side = midpoint_translate_side(g)
    for k in midpoint:
        if k % 2 == 0:
            midpoint[k] -= translate[k]
        else:
            a, b, c = side[k]
            side[k] = (a, b, -c)
    vertex = staircase_ngon_vertices(g, labels='zigzag')
    next_edge = dict((k, zfc((cfz(k) + 1) % n)) for k in range(1, 2*g+2))
    ref = reflect(g)
    exit = lambda d, p, i: hit_from_point_vector_to_line(d, p, side[i], V=V)
    S = substitution_words(g)
    w = S(list(map(Integer, w.split(',') if ',' in w else list(w))))
    if w not in p:
        augment(p, w, g=g)
    s = fold(p[w], g=g)  # the billiard sequence
    if print_length:
        print("%s: %s" %(w, len(s)))
    if print_cutting_word:
        print("%s: %s" %(w, s))
    v = vector_of_word(w, g=g)
    sp = midpoint[s[0]]
    if distance_along is not None:
        t = QQ(distance_along)
        starting_pt = t*vertex[s[0]] + (1-t)*vertex[next_edge[s[0]]]
        sp = starting_pt
    seg = []
    for i in s[1:]:
        ep = exit(sp, v, i)
        seg.append((sp, ep))
        sp = ep
        v = ref[i] * v
    # last step: close up path to very first start point
    ep = midpoint[s[0]]
    seg.append((sp, ep))
    res = SegmentSequence(segments=seg,
            rotation_order=None if return_rotation_order is None else len(s) // len(p[w]),
            cutting_word=None if return_cutting_word is None else s)
    return res


def matrix_L_to_P(g, base_ring=None):
    r"""
    """
    n = 2 * g + 1
    if base_ring is None:
        base_ring = AA
    if base_ring == RDF:
        pi_over_n = RDF.pi()/n
        c = pi_over_n.cos()
        s = pi_over_n.sin()
    else:
        pi_over_n = pi/n
        c = base_ring(cos(pi_over_n))
        s = base_ring(sin(pi_over_n))
    return matrix(base_ring, 2, [1, c, 0, s])


def draw_ngon_billiard_orbit(p, w, g='genus', v=None,
                             line_thickness=1,
                             adapt_line_thickness=True,
                             ngon_thickness=1,
                             print_cutting_word=False,
                             print_length=False,
                             border=True,
                             border_color='black',
                             line_color='black',
                             for_laser_cutting=False,
                             offset=(0, 0),
                             path_range=(0, 1),
                             distance_along=None,
                             start_from_edge=None,
                             evens=True,
                             odds=True,
                             kill_rotational_symmetry=False):
    """
    Return picture of regular pentagon with billiard orbit plotted.

    The picture is returned as a Sage graphics object.
    Aspect ratio is 1, no axes.

    The optional parameter ``offset`` is not taken into account if
    ``for_laser_cutting`` is `True`.

    The optional argument ``path_range`` lets you decide what range
    of the line segments of the orbit to actually plot. Useful for plotting
    an incomplete periodic trajectory (cheap way to pretend we're plotting
    a non-periodic trajectory).

    The optional parameter ``kill_rotational_symmetry`` can be ``False`` (default)
    or 1 or 2. Trajectories with rotational symmetry go through each midpoint twice,
    and 1 or 2 will make a different choice of how to cut.

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: G = Graphics()
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: draw_ngon_billiard_orbit(p, S([2, 1, 0]), g=g)
    """
    # TODO: fix the edge_reorder business
    n = 2*g + 1
    cfz = lambda k: cyclic_from_zigzag(g, k)
    edge_reorder = dict((k, dict((j, (cfz(j) - cfz(k)) % n)
                         for j in range(1, n + 1)))
                        for k in range(1, n + 1))
    billiard_segments = ngon_billiard_segments(p, w, g=g, v=v,
            print_cutting_word=print_cutting_word,
            print_length=print_length,
            distance_along=distance_along,
            return_rotation_order=True,
            return_cutting_word=True)
    seg = billiard_segments.segments
    billiard_word = billiard_segments.cutting_word
    rotation_order = billiard_segments.rotation_order
    # TODO: fix the kill_rotational_symmetry business
    if kill_rotational_symmetry == 1 and is_five_fold_symmetric_word(w):
        seg = seg[:len(seg)/5]
    if kill_rotational_symmetry == 2 and is_five_fold_symmetric_word(w):
        n = len(seg)/10
        seg = seg[n:3*n]
    if print_cutting_word:
        print(billiard_word)
    start_edge = billiard_word[0]
    aa = number_sequence(g, all_the_way=True)
    V = seg[0][0].parent()
    mm = midpoint_reverse_zigzag(g)
    ngon = [mm[0] - V((1, 0))]
    for m in mm[:1] + mm[1::2] + mm[:0:-2]:
        ngon.append(2*m - ngon[-1])
    center = mean(ngon[:-1])
    pi_over_n = RDF.pi()/n
    M = matrix(RDF, 2, [1, pi_over_n.cos(), 0, pi_over_n.sin()])
    ngon = [M * RDF_vector_from_nf_vector(x - center) for x in ngon]
    GG = Graphics()
    offset = vector(RDF, offset)
    if for_laser_cutting:
        line_thickness = 0.1
        line_color = 'black'
        GG += line2d([1.03 * p + offset for p in ngon], color='red',
                     aspect_ratio=1, axes=False, thickness=.1)
    if start_from_edge is not None and start_edge != start_from_edge:
        k = edge_reorder[start_edge][start_from_edge]
        angle = 2*k*pi_over_n
        M = matrix(RDF, 2, (cos(angle), -sin(angle), sin(angle), cos(angle)))*M
    if path_range != (0, 1):
        a, b = path_range
        seg = seg[floor(a * len(seg)):floor(b * len(seg))]
    f = RDF_vector_from_nf_vector
    seg = [(M * f(x - center), M * f(y - center)) for (x, y) in seg]
    if adapt_line_thickness:
        orbit_length = sum(norm(y - x) for x, y in seg)/norm(ngon[0])
        thickness_from_length = lambda x: 2/(1 + x/100)
        line_thickness = thickness_from_length(orbit_length)
    for k, (x, y) in enumerate(seg):
        if evens and k % 2 == 0:
            GG += line2d([x + offset, y + offset],
                 thickness=line_thickness, color=line_color)
        if odds and k % 2 == 1:
            GG += line2d([x + offset, y + offset],
                         thickness=line_thickness, color=line_color)
    if border:
        GG += line2d([p + offset for p in ngon], color=border_color,
                     aspect_ratio=1, axes=False, thickness=ngon_thickness)
        return GG


def compute_tiles(p, w, g='genus', v=None, show_progress=True, start_from_edge=4, split_by_parity=False):
    """
    Return tiles, i.e., complement of thickened periodic orbit in thickened ngon.

    TODO:

    - make it plot only polygons, not vertices (which show up as dots)
    --- to do this, do sum(tile.plot(hue=random(), point=False) for tile in tiles)

    Return picture of regular pentagon with billiard orbit plotted.

    The picture is returned as a string in svg syntax,
    describing the thickened periodic path.

    The optional argument ``path_range`` lets you decide what range
    of the line segments of the orbit to actually plot. Useful for plotting
    an incomplete periodic trajectory (cheap way to pretend we're plotting
    a non-periodic trajectory).

    The option split_by_parity can be used to get a two-coloring of the regions.

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: compute_tiles(p, S([1, 2, 0]), g=g)

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: tpos,tneg = compute_tiles(p, S([1, 0]), g=g, split_by_parity=True)
        sage: options = {'point': False, 'line': False, 'axes': False}
        sage: picpos = sum((plot(t, color='maroon', **options) for t in tpos), Graphics())
        sage: picneg = sum((plot(t, color='antiquewhite',  **options) for t in tneg), Graphics())
        sage: pic = picpos + picneg
        sage: pic.show()
    """
    edge_reorder = {
        # {current_start_edge: {desired_start_edge: how_much_to_rotate}}
        '1': {1: 0, 2: 4, 3: 1, 4: 3, 5: 2},
        '3': {1: 4, 2: 3, 3: 0, 4: 2, 5: 1},
        '5': {1: 3, 2: 2, 3: 4, 4: 1, 5: 0},
        '4': {1: 2, 2: 1, 3: 3, 4: 0, 5: 4},
        '2': {1: 1, 2: 0, 3: 2, 4: 4, 5: 3},
        }
    seg, billiard_word = compute_ngon_billiard_segments(p, w, g=g, v=v,
            return_billiard_word=True)
    start_edge = billiard_word[0]
    c, s = AA(cos(pi/5)), AA(sin(pi/5))
    M = matrix(AA, 2, [1, c, 0, s])
    V = VectorSpace(AA, 2)
    pentagon = [ (0, 0), (1, 0), (0, _u), (-_u, 1 + _u), (-_u, _u), (0, 0), ]
    pentagon = [ V(x) for x in pentagon ]
    center = 1/5 * sum(pentagon)
    pentagon = [M * (x - center) for x in pentagon] # now centered at zero
    polypos = [Polyhedron(pentagon)]
    polyneg = []
    if start_from_edge is not None and start_edge != start_from_edge:
        k = edge_reorder[start_edge][start_from_edge]
        angle = 2*k*pi/5
        c, s = AA(cos(angle)), AA(sin(angle))
        M = matrix(AA, 2, (c, -s, s, c))*M
    if show_progress:
        print('Number of segments is {}'.format(len(seg)))
    for k in range(len(seg)):
        if show_progress:
            print(k)
        x, y = seg[k]
        l = Polyhedron([M * (V(x) - center), M * (V(y) - center)])
        ppos = []
        pneg = []
        le = l.equations()[0]
        A = le.A()
        b = le.b()
        hl = Polyhedron(ieqs=[[b] + list(A)])
        hr = Polyhedron(ieqs=[[-b] + list(-A)])
        for p in polypos:
            pl = p.intersection(hl)
            pr = p.intersection(hr)
            if pl.n_vertices() > 2:
                ppos.append(pl)
            if pr.n_vertices() > 2:
                pneg.append(pr)
        for p in polyneg:
            pl = p.intersection(hl)
            pr = p.intersection(hr)
            if pl.n_vertices() > 2:
                pneg.append(pl)
            if pr.n_vertices() > 2:
                ppos.append(pr)
        polypos = ppos[:]
        polyneg = pneg[:]
    if show_progress:
        print()
        print('Number of polygons: {}'.format(len(polypos) + len(polyneg)))
    if split_by_parity:
          return polypos, polyneg
    return polypos + polyneg


def draw_thick_ngon_path(p, w, g='genus', v=None, line_thickness=1/20,
                             ngon_thickness=1/20, start_from_edge=4, show_progress=True):
    r"""
    Return billiard path plus outer pentagon, as an outline in svg.

    The picture is returned as a string in svg syntax,
    describing the thickened periodic path.

    .. TODO:

        - make sure the line thickness is in human-understandable units
        - from A,b do A,b+|A|*w/2
        - from Ax + b >= 0 do Ax + b + d/|A|

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: draw_thick_ngon_path(p, S([1, 2, 0]), g=g)

    """
    c, s = AA(cos(pi/5)), AA(sin(pi/5))
    a = 2*c
    M = matrix(AA, 2, [1, c, 0, s])
    V = VectorSpace(AA, 2)
    ngon = [ (0, 0), (1, 0), (0, a), (-a, 1 + a), (-a, a), (0, 0), ]
    ngon = [ V(x) for x in ngon ]
    center = 1/5 * sum(ngon)
    ngon = [M * (x - center) for x in ngon] # now centered at zero
    ngon = Polyhedron(ngon)
    polytopes = compute_tiles(p, w, g=g, v=v, show_progress=show_progress)
    thin_polytopes = []
    for k, p in enumerate(polytopes):
        if show_progress and k % 10 == 0:
            print('thinning polygon #{}'.format(k))
        q = Polyhedron(ieqs=[[i.b()-i.A().norm()*line_thickness/2]+list(i.A()) for i in p.inequalities()])
        if q.n_vertices() > 2:
           thin_polytopes.append(q)
    thick_ngon = Polyhedron(ieqs=[[i.b()+i.A().norm()*line_thickness/2]+list(i.A()) for i in ngon.inequalities()])
    # thin_polytopes.append(thick_ngon)
    return thin_polytopes, thick_ngon


def compute_necklace_orbit_segments(p, w, g='genus', v=None, verbose=False):
    """
    Orbit segments in slanted double polygon for substitution word ``w``

    INPUT:

    - ``p`` -- a dictionary mapping contraction words to cutting sequences

    - ``w`` -- a contraction word

    - ``v`` (default: `None`) -- ...

    OUTPUT:

    - A list of orbit segments (as pairs of points)

    EXAMPLE::

        sage: ...

    """
    # labels: '12' means side 1 of "double-pentagon" 2
    # pentagons are labeled 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
    # counterclockwise starting at the very bottom
    midpoint = { # the midpoints of the sides the trajectory leaves through
        '11': _V((0, -1/2)),
        '12': _V((0, _u + 1/2)),
        '13': _V((-2*_u - 1, 3*_u + 3/2)),
        '14': _V((-3*_u-2, 3*_u + 3/2)),
        '15': _V((-2*_u - 1, _u + 1/2)),
        '21': _V((-_u/2, _u + 1/2)),
        '22': _V((-3*_u/2 - 1, 3*_u + 3/2)),
        '23': _V((-7*_u/2 - 2, 4*_u + 5/2)),
        '24': _V((-7*_u/2 - 2, 3*_u + 3/2)),
        '25': _V((-3*_u/2 - 1, _u + 1/2)),
        '31': _V((-_u/2, _u/2)),
        '32': _V((-_u/2, 3/2*_u+1)),
        '33': _V((-5*_u/2-1, 7*_u/2+2)),
        '34': _V((-7*_u/2-2, 7*_u/2+2)),
        '35': _V((-5*_u/2-1, 3/2*_u+1)),
        '41': _V((1/2, _u/2)),
        '42': _V((-_u -1/2, 5*_u/2 + 1)),
        '43': _V((-3*_u-3/2, 7*_u/2+2)),
        '44': _V((-3*_u-3/2, 5*_u/2+1)),
        '45': _V((-_u-1/2, _u/2)),
        '51': _V((-_u-1/2, _u)),
        '52': _V((-_u-1/2, 2*_u + 1)),
        '53': _V((-3*_u-3/2, 4*_u+2)), # changed -3*_u/2 to -5*...
        '54': _V((-4*_u-5/2, 4*_u+2)),
        '55': _V((-3*_u-3/2, 2*_u+1)),
        }
    translate = { # translation vectors from the leaving edge to the entering edge -- all correct
        '11': _V((-4*_u-2, 3*_u+2)),
        '12': _V((-4*_u-2, 3*_u+2)),
        '13': _V((0, 0)),
        '14': _V((2*_u+2, -2*_u-1)),
        '15': _V((0, 0)),
        '21': _V((0, 0)),
        '22': _V((_u+1, -3*_u-2)),
        '23': _V((_u+1, -3*_u-2)),
        '24': _V((0, 0)),
        '25': _V((-_u, 2*_u+1)),
        '31': _V((0, 0)),
        '32': _V((-_u-1, -_u-1)),
        '33': _V((-_u-1, -_u-1)),
        '34': _V((0, 0)),
        '35': _V((_u, _u)),
        '41': _V((-3*_u-2, _u + 1)),
        '42': _V((-3*_u-2, _u + 1)),
        '43': _V((0, 0)),
        '44': _V((2*_u+1, -_u)),
        '45': _V((0, 0)),
        '51': _V((-2*_u-1, 2*_u+2)),
        '52': _V((0, 0)),
        '53': _V((3*_u+2, -4*_u-2)),
        '54': _V((3*_u+2, -4*_u-2)),
        '55': _V((0, 0)),
        }
    side = { # the 3 coordinates are a,b,c in ax+by=c, the line supporting the edge that is hit
        '11': (1, 0, 0),
        '12': (1, 0, 0),
        '13': (1, 0, -2*_u-1),
        '14': (1, 0, -3*_u-2),
        '15': (1, 0, -2*_u-1),
        '21': (1, _u, _u+1),
        '22': (1, _u, 3*_u+2),
        '23': (1, _u, 3*_u+2),
        '24': (1, _u, _u+1),
        '25': (1, _u, 0),
        '31': (1, 1, 0),
        '32': (1, 1, _u+1),
        '33': (1, 1, _u+1),
        '34': (1, 1, 0),
        '35': (1, 1, -_u),
        '41': (_u, 1, _u),
        '42': (_u, 1, _u),
        '43': (_u, 1, -_u-1),
        '44': (_u, 1, -2*_u-2),
        '45': (_u, 1, -_u-1),
        '51': (0, 1, _u),
        '52': (0, 1, 2*_u+1),
        '53': (0, 1, 4*_u+2),
        '54': (0, 1, 4*_u+2),
        '55': (0, 1, 2*_u+1),
        }
    dpcopy = { # which double pentagon copy we go to when exiting that side
        '11': '4',
        '12': '3',
        '13': '2',
        '14': '1',
        '15': '5',
        '21': '2',
        '22': '1',
        '23': '5',
        '24': '4',
        '25': '3',
        '31': '1',
        '32': '5',
        '33': '4',
        '34': '3',
        '35': '2',
        '41': '5',
        '42': '4',
        '43': '3',
        '44': '2',
        '45': '1',
        '51': '3',
        '52': '2',
        '53': '1',
        '54': '5',
        '55': '4',
    }
    cheatsheet = { # to enter side i of dp 1, need to exit from which dp
        '1': '4',
        '2': '2',
        '3': '1',
        '4': '5',
        '5': '3',
    }
    def endpoint(sp, v, s):
        # sp: starting point, v: vector, s: exit side
        # since we know the cutting sequence in advance, a starting point, a direction vector,
        # and the exit side label, which records which side of which pentagon,
        # and thus the equation of the line we will hit, are enough
        # to compute the exit point.
        a, b = sp
        c, d = v
        aa, bb, xx, cc, dd, yy = side[s] + (-d, c, b * c - a * d)
        delta = (aa * dd - bb * cc)
        return _V(((xx * dd - bb * yy) / delta, (aa * yy - xx * cc) / delta))
    if w not in p:
        augment(p, w)
    s = p[w] # the double pentagon cutting sequence
    if verbose:
        print("%s: %s" %(w, p[w]))
    v = vector_of_word(w, v)
    if has_five_fold_symmetry_vec(v):
        s *= 5
    # want to start in dpcopy 1 entering from side s[0]
    pent = '1'
    j = cheatsheet[s[0]]
    # final end point
    fep = midpoint[s[0] + j]
    # initial starting point
    sp = fep + translate[s[0] + j]
    seg = []
    if verbose:
        print('s = {}'.format(s))
        print('v = {}'.format(v))
    for i in s[1:]:
        ep = endpoint(sp, v, i + pent)
        if verbose:
            print('i = {}'.format(i))
            print('pent = {}'.format(pent))
            print('sp = {}'.format(sp))
            print('ep = {}'.format(ep))
        seg.append((sp, ep))
        sp = ep + translate[i + pent]
        pent = dpcopy[i + pent]
    seg.append((sp, fep))
    return seg


def draw_necklace_orbit(p, w, g='genus', v=None,
                        line_thickness=1,
                        polygon_thickness=1,
                        verbose=False):
    """
    Return picture of double regular polygon with orbit plotted

    The picture is returned as a Sage graphics object.
    Double polygon in red, line segments for the orbit in blue.
    Aspect ratio is 1, no axes.

    EXAMPLE::

        sage: # load('gorbits.sage')
        sage: g = 2
        sage: n = 2*g + 1
        sage: S = substitution_words(g)
        sage: C = cutting_words(g)
        sage: p = p12 = {S(): C([1, 2])}
        sage: draw_dp_orbit(p, S([1, 2, 0]), g=g)
    """
    seg = compute_necklace_orbit_segments(p, w, v=v, verbose=verbose)
    pi5 = RDF.pi()/5
    M = matrix(RDF, 2, [1, pi5.cos(), 0, pi5.sin()])
    necklace_out = [
                (0,0), (1,0), (0,_u), (0, 1+_u), (-_u, 1+2*_u),
                (-1-_u,1+3*_u), (-1-2*_u,2+3*_u), (-1-3*_u, 2+4*_u),
                (-2-3*_u, 2+4*_u), (-2-4*_u, 3+4*_u),
                (-2-4*_u,2+4*_u), (-3-4*_u,2+4*_u), (-2-4*_u,2+3*_u),
                (-2-4*_u, 1+3*_u), (-2-3*_u,1+2*_u), (-1-3*_u,1+_u),
                (-1-2*_u,_u), (-1-_u,0), (-_u, 0), (0, -1), (0, 0)]
    necklace_in = [
            (-_u, _u), (-_u, _u + 1), (-_u - 1, 2*_u + 1),
            (-2*_u - 1, 3*_u + 1), (-3*_u - 1, 3*_u + 2),
            (-3*_u - 2, 3*_u + 2), (-3*_u - 2, 3*_u + 1),
            (-3*_u - 1, 2*_u + 1), (-2*_u - 1, _u + 1),
            (-_u - 1, _u), (-_u, _u)]
    necklace_out = map(lambda x: M * x, [ _V(x) for x in necklace_out ])
    necklace_in = map(lambda x: M * x, [ _V(x) for x in necklace_in ])
    GG = Graphics()
    i = 0
    for x, y in seg:
        GG += line2d([M * x, M * y], thickness=line_thickness)
        if verbose:
            GG += text(str(i) + '>', M * (9/10*x + 1/10*y))
            GG += text('>' + str(i), M * (1/10*x + 9/10*y))
        i += 1
    GG += line2d(necklace_out, color='red', aspect_ratio=1, axes=False,
                 thickness=polygon_thickness)
    GG += line2d(necklace_in, color='red', aspect_ratio=1, axes=False,
                 thickness=polygon_thickness)
    for i in range(10):
        GG += line2d([necklace_in[i], necklace_out[2*i]],
                    color='red', aspect_ratio=1, axes=False,
                    thickness=polygon_thickness)
    return GG


def draw_2g_orbits(p, w, g=4, v=None, **options):
    S = substitution_words(g)
    if w not in S:
        w = S(list(map(Integer, w.split(',') if ',' in w else list(w))))
    pic = draw_ngon_billiard_orbit(p, w + S([0]), g=g, v=v)
    xyrange = pic.get_axes_range()
    xoffset = 1.125 * (xyrange['xmax'] - xyrange['xmin'])
    yoffset = 1.125 * (xyrange['ymax'] - xyrange['ymin'])
    ngons = [pic]
    for x in srange(1, 2*g):
        ww = w + S([x])
        args = (p, w + S([x]))
        a, b = x.quo_rem(4)
        opts = {'g': g, 'v': v, 'offset': (b * xoffset, -a * yoffset)}
        opts.update(options)
        ngons.append(draw_ngon_billiard_orbit(*args, **opts))
    return sum(ngons, Graphics())


def four_orbits_alt(p, w, v=None, line_thickness=1, pentagon_thickness=1,
                print_cutting_word=False, print_lengths=False, save=False):
    if save:
        graphics_array(
            [draw_pentagon_billiard_orbit(p, w+x, v=v,
                                          line_thickness=line_thickness,
                                          pentagon_thickness=pentagon_thickness,
                                          print_cutting_word=print_cutting_word,
                                          print_length=print_lengths)
             for x in '0123']).save('pentagons-{}-{}.png'.format(save, w), axes=False)
    a, b, c, d = [draw_pentagon_billiard_orbit(p, w+x, v=v,
                                      line_thickness=line_thickness,
                                      pentagon_thickness=pentagon_thickness,
                                      print_cutting_word=print_cutting_word,
                                      print_length=print_lengths, offset=(1.8*int(x), 0))
         for x in '0123']
    return a + b + c + d


def eight_orbits(p, pp, w, v=None, vv=None, line_thickness=0.2, pentagon_thickness=0.8):
    return graphics_array(
        [draw_pentagon_billiard_orbit(p,w+x,v=v,line_thickness=line_thickness,
                                      pentagon_thickness=pentagon_thickness)
         for x in '0123']
        + [draw_pentagon_billiard_orbit(pp,w+x,v=vv,line_thickness=line_thickness,
                                        pentagon_thickness=pentagon_thickness)
           for x in '0123'],
        2, 4).show(axes=False)


def eight_orbits_alt(p, pp, w, v=None, vv=None, line_thickness=0.2, pentagon_thickness=0.8):
        a, b, c, d = [draw_pentagon_billiard_orbit(p,w+x,v=v,line_thickness=line_thickness,
                                      pentagon_thickness=pentagon_thickness, offset=(1.8*int(x), 0))
         for x in '0123']
        aa, bb, cc, dd = [draw_pentagon_billiard_orbit(pp,w+x,v=vv,line_thickness=line_thickness,
                                        pentagon_thickness=pentagon_thickness, offset=(1.8*int(x), -1.7))
           for x in '0123']
        return a + b + c + d + aa + bb + cc + dd


def has_five_fold_symmetry_vec(v):
    a, b = v[0]
    c, d = v[1]
    return bool(((d - b) + 2*(c - a)) % 5)


def is_five_fold_symmetric_word(w):
    M = MatrixSpace(ZZ, 4)
    a = M((1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1))
    b = M((0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1))
    c = M((0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1))
    d = M((1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1))
    mm = [a, b, c, d]
    onezero = vector(ZZ, (1, 0, 0, 0))
    a, b, c, d = prod(mm[int(i)] for i in w[::-1])*onezero
    return bool(((d - b) + 2*(c - a)) % 5)
