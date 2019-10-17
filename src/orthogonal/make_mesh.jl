
rdata = joinpath(rdata,"sint_nero")
OUT = [joinpath(rdata,"4_mesh.tsv"),joinpath(rdata,"5_geom.tsv"),joinpath(rdata,"6_wellCon.tsv")]
xy = collect(Iterators.product(10:20:1000, 10:20:1000))
xy = convert(Array{Tuple{Float64,Float64},2},xy)
xy = map(x->[x[1],x[2]],xy)[:]
xy = xy/1000
wxy = [[250 250],
       [750 750]]

make_mesh(xy,OUT[1])
make_geom(xy,OUT[2])
make_wellCon(wxy,xy,OUT[3])

function make_wellCon(wxy,xy,OUT)
    nw = length(wxy)
    w1 = zeros(Int64,nw)
    for i=1:nw
        w1[i] = findmin(map(x->sum((x.-wxy[i]).^2),xy[:]))[2]
    end

    w2 = collect(1:nw)
    ioW = Base.open(OUT,"w");
    writeToPipe(ioW, w2, w1)
    close(ioW)
end

function make_geom(xy,OUT)
    n = length(xy)
    id = collect(1:n)
    id50 = reshape(id,50,50)
    r = map(x->[x-1,x+1,x-50,x+50],id50)
    r = map(x->x[0 .<x.<=2500],r)[:]
    c = id50[:]
    b = map(x->fill(20.,length(x)),r)
    area_edge = map(x->fill(100.,length(x)),r)
    l = map(x->fill(20.,length(x)),r)
    ztop = fill(1000.,n)
    zbot = fill(1010.,n)
    area = fill(400.,n)

    r = Ar2Str(r)
    b = Ar2Str(b)
    area_edge = Ar2Str(area_edge)
    l = Ar2Str(l)

    ioW = Base.open(OUT,"w");
    writeToPipe(ioW, c, r, b, area_edge, l, ztop, zbot, area)
    close(ioW)
end

function Ar2Str(r)
    r = map(x->"$x",r)
    r= map(x->replace(x, "[" => "\""),r)
    r = map(x->replace(x, "]" => "\""),r)
    return r
end

function make_mesh(xy,OUT)
    n = length(xy)
    id = collect(1:length(xy))
    xy =  convert(Array{Tuple{Float64,Float64},2},xy)
    bnd = zeros(Int32,size(xy));
    id50 = reshape(id,50,50)
    bnd[1,1:end-1] .= id50[1,2:end];
    bnd[1:end-1,end] .= id50[2:end,end];
    bnd[end,2:end] .= id50[end,1:end-1];
    bnd[2:end,1] .= id50[1:end-1,1];

    X = map(x->x[1],xy[:]);
    Y = map(x->x[2],xy[:]);
    bnd = bnd[:]

    vxB = Matrix{String}(undef,size(xy))
    vyB = Matrix{String}(undef,size(xy))
    vxB = map(x->"$(x[1]-10), $(x[1]+10)",xy)
    vyB = map(x->"$(x[2]-10), $(x[2]+10)",xy)

    vxB[1] = "0, $(vxB[1])"
    vxB[50] = "$(vxB[50]), 1000"
    vxB[2451] = "0, $(vxB[2451])"
    vxB[2500] = "$(vxB[2500]), 1000"

    vyB[1] = "0, $(vyB[1])"
    vyB[50] = "0, $(vyB[50])"
    vyB[2451] = "$(vyB[2451]), 1000"
    vyB[2500] = "$(vyB[2500]), 1000"

    for i=1:2500
        vxB[i] = "\"$(vxB[i])\""
        vyB[i] = "\"$(vyB[i])\""
    end

    vxB[2:end-1,2:end-1].="\\N"
    vyB[2:end-1,2:end-1].="\\N"

    hz = fill("1",n)
    ioW1 = Base.open(OUT,"w");
    writeToPipe(ioW1, id, X, Y, bnd, hz, vxB[:], vyB[:])
    close(ioW1)
end
