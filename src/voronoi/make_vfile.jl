function make_vfile(OUT,wi,wxy,rc,Cv, exy)
    xyc_ab, xy_ab = getFromCV(Cv);
    ic_ab, xyc = get_index(xyc_ab); #индексы центров и координаты
    ie_ab, xye = get_index(xy_ab);  #индексы границ и координаты

    i_bnd = findall((x->all(x==(0.,0.))).(xyc))[1]; #Индекс отвечающий за границу
    i_ed = (x->x[1]).(findall(ic_ab.==i_bnd)); #Индексы связей, которые являются границей;

    bnd_ind = zeros(Int32,length(i_ed),2)
    ebo = eage_order(xye,ie_ab,i_bnd,ic_ab); #порядок индексов границы иджей из xye
    pbo = point_opder(ebo,ie_ab,ic_ab,i_bnd); #порядок индексов границы центров из xyс

    ib = all(indexin(ie_ab,ebo).!=nothing,dims=2)[:];
    ib1 = ic_ab[ib,1]
    ib2 = ie_ab[ib,:]

    xyc = xyc[setdiff(1:length(xyc),i_bnd)]

    xyc = mk_decmprs(xyc,exy)
    xye = mk_decmprs(xye,exy)

    make_mesh(xyc,xye,pbo,ib1,ib2,OUT[1])
    make_geom(xyc,rc,ie_ab, ic_ab, xye, OUT[2])
    make_wellCon(wi,wxy,xyc,OUT[3])

    return true
end

function eage_order(xye,ie_ab,i_bnd,ic_ab)
    next_p = zeros(Int32,length(xye))

    #индекс первой точки обхода
    #next_p[1] = findfirst((x->all(x.==(1.25,1.25))).(xye));
    next_p[1] = findmin(map(x->sum(x.^2),xye))[2]

    # ia = (x->x[1]).(findall(ie_ab.==next_p[1]));
    # next_p[2] = setdiff(ie_ab[ia[1],:],ind_1)[1]

    fl = true
    k=0
    while fl
        k+=1;
        #println(xye[next_p[k]])
        ia = (x->x[1]).(findall(ie_ab.==next_p[k])); #индексы связей где встречается эта точка
        #println("$k $(length(ia))")
        ia=ia[any(ic_ab[ia,:].==i_bnd,dims=2)[:]]

        v = setdiff(ie_ab[ia,:],next_p[1:k]);
        if length(v)>0
            next_p[k+1] = v[1]
        else
            fl = false
        end

    end
    return next_p[1:k]
end

function point_opder(ebo,ie_ab,ic_ab,i_bnd)
    pbo = zeros(Int32,length(ebo))

    ia = findall(ebo[1].==ie_ab)
    ia = (x->x[1]).(ia)
    ib = ic_ab[ia,:]
    ti = setdiff(unique(ib),i_bnd);
    pbo[1] = ti[1]

    k=1;
    for v in Iterators.rest(ebo,2)
        ia = findall(v.==ie_ab)
        ia = (x->x[1]).(ia)
        ib = ic_ab[ia,:]
        ti = setdiff(unique(ib),i_bnd);
        g = setdiff(ti,pbo[k])
        if length(g)>0
            k+=1;
            pbo[k] = g[1]
        end
    end
    return pbo[1:k]
end

function get_index(xy_ab)
    xyab = collect(zip(view(xy_ab,:,1:2)[:],view(xy_ab,:,3:4)[:]))
    xy = unique(xyab); #xy =xy[.!(x->all(x==(0.,0.))).(xy)]
    ia = indexin(xyab,xy)
    return reshape(ia,size(xy_ab,1),2), xy
end

function make_mesh(xy,xye,pbo,ib1,ib2,OUT)
    n = length(xy)
    id = collect(1:length(xy))
    #xy =  convert(Array{Tuple{Float64,Float64},2},xy)
    bnd = zeros(Int32,size(xy));
    bnd[pbo[1:end-1]].=pbo[2:end]

    X = map(x->x[1],xy[:]);
    Y = map(x->x[2],xy[:]);
    bnd = bnd[:]

    vxB = Vector{String}(undef,length(xy))
    vyB = Vector{String}(undef,length(xy))

    vxB[:].="\\N"
    vyB[:].="\\N"

    #ia = indexin(ib1,pbo)
    for i in pbo
        ia = findall(i.==ib1)
        txy = xye[o1(ib2[ia,:])]
        vxB[i] =  join(((x->"$(x[1])").(txy)),", ")
        vyB[i] =  join(((x->"$(x[2])").(txy)),", ")
    end

    hz = fill("1",n)
    ioW1 = Base.open(OUT,"w");
    writeToPipe(ioW1, id, X, Y, bnd, hz, vxB[:], vyB[:])
    close(ioW1)
end

function make_wellCon(wi,wxy,xy,OUT)
    nw = length(wi)
    w1 = zeros(Int64,nw)
    wxy = mat2vec(wxy)
    for i=1:nw
        w1[i] = findmin(map(x->sum((x.-wxy[i]).^2),xy[:]))[2]
    end
    w2 = wi;#collect(1:nw)
    w21 = unique(hcat(w2,w1),dims=1);
    ioW = Base.open(OUT,"w");
    writeToPipe(ioW, view(w21,:,1), view(w21,:,2))
    close(ioW)
end

function make_geom(xy,rc,ie_ab, ic_ab, xye,OUT)
    n = length(xy)
    id = collect(1:n)

    #id50 = reshape(id,50,50)
    brc, lrc = distFromRC(rc,ie_ab, ic_ab, xy,xye);
    c = id[:]
    r = Vector(undef,n)
    b = Vector(undef,n)
    l = Vector(undef,n)
    area_edge = Vector(undef,n)
    area = zeros(Float64,n)
    for ic in c
        ni = (x->x[1]).(findall(rc.==ic))
        r[ic] = circshift(rc,(0,1))[rc.==ic]
        b[ic] = brc[ni]
        l[ic] = lrc[ni]
        eo = o2(ie_ab[any(ic.==ic_ab,dims=2)[:],:])
        area[ic] = polyaria(vcat(map(x->[x[1],x[2]], xye[eo])'...))
        area_edge[ic] = b[ic]*10;
    end

    ztop = fill(1000.,n)
    zbot = fill(1010.,n)


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
function dist(vxy)
    sqrt(sum((vxy[1].-vxy[2]).^2))
end

function distFromRC(rc,ie_ab, ic_ab, xyc, xye)
    db = zeros(Float64,size(rc,1))
    dl = zeros(Float64,size(rc,1))
    for i=1:size(rc,1)
        poe = ie_ab[findfirst(all((rc[i,1].==ic_ab) .| (rc[i,2].==ic_ab),dims=2)[:]),:]
        db[i]=dist(xye[poe])
        dl[i]=dist(xyc[rc[i,:]])
    end
    return db, dl
end
