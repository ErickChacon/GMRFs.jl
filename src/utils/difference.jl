"""
    difference(d, order = 1, circular = false)

Return a sparse `d×n` difference matrix D for a domain `d`, where `n` represents the
domain's size and `d` is the number of differences (restrictions or increments).

Each row in D has coefficients that define the restrictions of specified `order`  with
respect to the `n` elements. Supported domain types include `CartesianGrid` or
`SimpleGraph`, with an option for CartesianGrid to treat the grid as circular (`circular =
true`), allowing neighbors to wrap around the grid's edges like a toroidal topology.
"""
function difference(g::CartesianGrid; order = 1, circular = false)
    Dim = paramdim(g)
    if Dim == 1
        return difference1(g; order = order, circular = circular)
    elseif Dim == 2
        return difference2(g; order = order, circular = circular)
    else
        throw(ErrorException("not implemented for Dim > 2"))
    end
end

function difference1(g::CartesianGrid; order = 1, circular = false)
    n = nelements(g)

    if order == 1 && !circular
        return spdiagm(n-1, n, 0 => fill(-1, n-1), 1 => fill(1, n-1))
    end

    if order == 1 && circular
        return spdiagm(0 => fill(-1, n), 1 => fill(1, n-1), -(n-1) => fill(1, 1))
    end

    if order == 2 && !circular
        return spdiagm(n-2, n,
                       0 => fill(1, n-2), 1 => fill(-2, n-2), 2 => fill(1, n-2))
    end

    if order == 2 && circular
        return spdiagm(0 => fill(1, n), 1 => fill(-2, n-1), 2 => fill(1, n-2),
                       -(n-2) => fill(1, 2), -(n-1) => fill(-2, 1))
    end

    throw(ErrorException("not implemented"))
end

function difference2(g::CartesianGrid; order = 1, circular = false)
    n1, n2 = size(g)
    n = nelements(g)

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        (j - 1) * n1 + i
    end

    if order == 1
        if circular
            # neighbors to the right (→)
            Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
            Jr = [ij_to_k(mod1(i + 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
            # neighbors to the top (↑)
            It = Ir
            Jt = [ij_to_k(i, mod1(j + 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        else
            # neighbors to the right (→)
            Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
            Jr = [ij_to_k(i + 1, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
            # neighbors to the top (↑)
            It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
            Jt = [ij_to_k(i, j + 1, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
        end
        return sparse(
            repeat(1:(length(Ir) + length(It)), 2),
            vcat(Ir, It, Jr, Jt),
            repeat([-1, 1], inner = length(Ir) + length(It))
           )
    end

    if order == 2 && circular
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the right (→)
        Jr = [ij_to_k(mod1(i + 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the top (↑)
        Jt = [ij_to_k(i, mod1(j + 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the left (←)
        Jl = [ij_to_k(mod1(i - 1, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the bottom (↓)
        Jb = [ij_to_k(i, mod1(j - 1, n2), n1, n2) for i = 1:n1 for j = 1:n2]
        return sparse(
            repeat(1:length(Ir), 5),
            vcat(Ir, Jr, Jt, Jl, Jb),
            repeat([-4, 1, 1, 1, 1], inner = length(Ir))
           )
    end

    # TODO: need improvement?
    if order == 2 && !circular
        # inside cells
        Ir = [ij_to_k(i, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        Jir = [ij_to_k(i + 1, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        Jit = [ij_to_k(i, j + 1, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        Jil = [ij_to_k(i - 1, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        Jib = [ij_to_k(i, j - 1, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]
        # sides cells
        Isl = [ij_to_k(1, j, n1, n2) for j = 2:(n2-1)]
        Jsl1 = [ij_to_k(1 + 1, j, n1, n2) for j = 2:(n2-1)]
        Jsl2 = [ij_to_k(1, j + 1, n1, n2) for j = 2:(n2-1)]
        Jsl3 = [ij_to_k(1, j - 1, n1, n2) for j = 2:(n2-1)]
        Isr = [ij_to_k(n1, j, n1, n2) for j = 2:(n2-1)]
        Jsr1 = [ij_to_k(n1 - 1, j, n1, n2) for j = 2:(n2-1)]
        Jsr2 = [ij_to_k(n1, j + 1, n1, n2) for j = 2:(n2-1)]
        Jsr3 = [ij_to_k(n1, j - 1, n1, n2) for j = 2:(n2-1)]
        Ist = [ij_to_k(i, n2, n1, n2) for i = 2:(n1-1)]
        Jst1 = [ij_to_k(i + 1, n2, n1, n2) for i = 2:(n1-1)]
        Jst2 = [ij_to_k(i - 1, n2, n1, n2) for i = 2:(n1-1)]
        Jst3 = [ij_to_k(i, n2 - 1, n1, n2) for i = 2:(n1-1)]
        Isb = [ij_to_k(i, 1, n1, n2) for i = 2:(n1-1)]
        Jsb1 = [ij_to_k(i + 1, 1, n1, n2) for i = 2:(n1-1)]
        Jsb2 = [ij_to_k(i - 1, 1, n1, n2) for i = 2:(n1-1)]
        Jsb3 = [ij_to_k(i, 1 + 1, n1, n2) for i = 2:(n1-1)]
        Is = vcat(Isl, Isr, Ist, Isb)
        Js = vcat(Jsl1, Jsr1, Jst1, Jsb1,
                  Jsl2, Jsr2, Jst2, Jsb2,
                  Jsl3, Jsr3, Jst3, Jsb3)
        # corners cells
        Ic1 = ij_to_k(1, 1, n1, n2)
        Jc1 = [ij_to_k(1 + 1, 1, n1, n2), ij_to_k(1, 1 + 1, n1, n2)]
        Ic2 = ij_to_k(1, n2, n1, n2)
        Jc2 = [ij_to_k(1 + 1, n2, n1, n2), ij_to_k(1, n2 - 1, n1, n2)]
        Ic3 = ij_to_k(n1, 1, n1, n2)
        Jc3 = [ij_to_k(n1 - 1, 1, n1, n2), ij_to_k(n1, 1 + 1, n1, n2)]
        Ic4 = ij_to_k(n1, n2, n1, n2)
        Jc4 = [ij_to_k(n1 - 1, n2, n1, n2), ij_to_k(n1, n2 - 1, n1, n2)]
        Ic = vcat(Ic1, Ic2, Ic3, Ic4)
        Jc = vcat(Jc1[1], Jc2[1], Jc3[1], Jc4[1], Jc1[2], Jc2[2], Jc3[2], Jc4[2])
        return sparse(
            vcat(repeat(1:length(Ir), 5),
                 repeat((1:length(Is)) .+ length(Ir), 4),
                 repeat((1:length(Ic)) .+ length(Ir) .+ length(Is), 3),
                ),
            vcat(Ir, Jir, Jit, Jil, Jib, Is, Js, Ic, Jc),
            vcat(repeat([-4, 1, 1, 1, 1], inner = length(Ir)),
                 repeat([-3, 1, 1, 1], inner = length(Is)),
                 repeat([-2, 1, 1], inner = length(Ic)),
                )
           )
    end

    throw(ErrorException("not implemented"))
end

function difference(g::SimpleGraph; order = 1)

    if order == 1
        n = ne(g)
        In = [e.src for e in edges(g)]
        Jn = [e.dst for e in edges(g)]
        return sparse(
            repeat(1:n, 2),
            vcat(In, Jn),
            repeat([-1, 1], inner = n)
        )
    elseif order == 2
        D = difference(g)
        return -D'D
    else
        throw(ErrorException("not implemented"))
    end

end


