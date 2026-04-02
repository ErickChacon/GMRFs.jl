"""
    structure(d; δ = 0, order = 1, circular = false)

Return a sparse `n×n` structure matrix S for a domain `d`, where `n` represents the
domain's size.

The structure matrix is proportional to the precision matrix (`Q = κS`) of a Gaussian
Markov Random Field (GMRF). It is also related to the difference matrix (`S = D'D`), which
is used to define a GMRF through conditional increments. Higher `order` will lead to a
smoother GMRF. An optional paramater `δ` can be defined such as `S = D'D + δ*I`, where D
represents the difference matrix. Supported domain types include `CartesianGrid` and
`SimpleGraph`, with an option for CartesianGrid to treat the grid as circular (`circular =
true`), allowing neighbors to wrap around the grid's edges like a toroidal topology.
"""
function structure(g::CartesianGrid; δ = 0, order = 1, circular = false)
    D₁ = difference(g, circular = circular)
    if order == 1
        S = D₁'D₁
    elseif order > 1
        S = D₁'D₁ * structure(g, order = order - 1, circular = circular)
    end
    S + δ * I
end

function structure(g::SimpleGraph; δ = 0, order = 1)
    D₁ = difference(g, order = 1)
    if order == 1
        S = D₁'D₁
    elseif order > 1
        S = D₁'D₁ * structure(g, order = order - 1)
    end
    S + δ * I
end

"""
    structure_base(d::CartesianGrid; δ = 0, order = 1)

Return a sparse base matrix S for a Cartesian grid `d`, assuming a circular neighboring
structure where neighbors wrap around the grid's edges in a toroidal fashion.

Due to the circular topology, the resulting structure matrix is circulant, obviating the
need for a full structure matrix computation, and instead focusing on the base of this
circulant matrix. Optional parameters `δ` and `order` are supported, offering the same
interpretation as in the `structure` matrix.
"""
function structure_base(g::CartesianGrid; δ = 0, order = 1)
    Dim = paramdim(g)

    if Dim == 1
        return structure_base1(g; δ = δ, order = order)
    elseif Dim == 2
        return structure_base2(g; δ = δ, order = order)
    else
        throw(ErrorException("not implemented for Dim > 2"))
    end
end

function structure_base1(g::CartesianGrid; δ = 0, order = 1)
    n = nelements(g)
    if order == 1
        base = sparsevec([1, 2, n], [2.0, -1, -1])
    elseif order > 1
        w = centered([-1, 2, -1])
        base = imfilter(structure_base(g, order = order - 1), w, "circular")
    end
    base[1,1] += δ
    base
end

function structure_base2(g::CartesianGrid; δ = 0, order = 1)
    if order == 1
        (n1, n2) = size(g)
        base = sparse([1, 1, 2, 1, n2], [1, 2, 1, n1, 1], [4.0, -1, -1, -1, -1])
    elseif order > 1
        w = centered(sparse([0 -1 0; -1 4.0 -1; 0 -1 0]))
        base = imfilter(structure_base(g, order = order - 1), w, "circular")
    end
    base[1,1] += δ
    base
end
