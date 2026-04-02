"""
    adjacency(d, circular = false)

Return a sparse `n×n` adjacency matrix A for a domain `d`, where `n` represents the
domain's size.

The value of A[i,j] is `true` if the elements `i` and `j` are first-order neighbors;
otherwise, it is false. Supported domain types include `CartesianGrid`, `GeometrySet`, or
`SimpleGraph`, with an option for CartesianGrid to treat the grid as circular (`circular =
true`), allowing neighbors to wrap around the grid's edges like a toroidal topology.
"""
function adjacency(d::CartesianGrid; circular::Bool = false)
    Dim = paramdim(d)
    if Dim == 1
        adjacency1(d; order = 1, circular = circular)
    elseif Dim == 2
        adjacency2(d; order = 1, circular = circular)
    else
        throw(ErrorException("adjacency is not implemented for Dim > 2"))
    end
end

function adjacency1(d::CartesianGrid; order::Int = 1, circular::Bool = false)
    n = nelements(d)

    # neighbors to the right (→)
    if circular
        Ir = 1:n
        Jr = mod1.(Ir .+ order, n)
    else
        Ir = 1:(n - order)
        Jr = (1 + order):n
    end

    # sparse adjacency matrix
    A = sparse(Ir, Jr, true, n, n)
    A .| A'
end

function adjacency2(d::CartesianGrid; order::Int = 1, circular::Bool = false)
    n1, n2 = size(d)
    n = nelements(d)

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        (j - 1) * n1 + i
    end

    if circular
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
        Jr = [ij_to_k(mod1(i + order, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the top (↑)
        It = Ir
        Jt = [ij_to_k(i, mod1(j + order, n2), n1, n2) for i = 1:n1 for j = 1:n2]
    else
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        Jr = [ij_to_k(i + order, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        # neighbors to the top (↑)
        It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
        Jt = [ij_to_k(i, j + order, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
    end

    # sparse adjacency matrix
    A = sparse(vcat(Ir, It), vcat(Jr, Jt), true, n, n)
    A .| A'
end

function adjacency(d::GeometrySet)
    # initialize
    n = length(d)
    rows = Int64[]
    cols = Int64[]

    # neighbors
    for i = 1:n
        for j = (i+1):n
            if intersects(d[i], d[j])
                push!(rows, i)
                push!(cols, j)
            end
        end
    end

    # sparse adjacency matrix
    A = sparse(rows, cols, true, n, n)
    A .| A'
end

adjacency(d::SimpleGraph) = adjacency_matrix(d, Bool)
