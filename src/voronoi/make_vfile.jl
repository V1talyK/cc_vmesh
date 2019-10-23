function make_vfile(OUT,rc,Cv)
    xyc_ab, xy_ab = getFromCV(Cv1);
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
    make_mesh(xyc,pbo,OUT[1])
    make_geom(xy,OUT[2])
    make_wellCon(wxy,xy,OUT[3])

    return
end

function eage_order(xye,ie_ab,i_bnd,ic_ab)
    next_p = zeros(Int32,length(xye))
    next_p[1] = findfirst((x->all(x.==(1.,1.))).(xye)); #индекс первой точки обхода

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

function make_mesh(xy,pbo,OUT)
    n = length(xy)
    id = collect(1:length(xy))
    #xy =  convert(Array{Tuple{Float64,Float64},2},xy)
    bnd = zeros(Int32,size(xy));
    bnd[pbo[1:end-1]].=pbo[2:end]

    X = map(x->x[1],xy[:]);
    Y = map(x->x[2],xy[:]);
    bnd = bnd[:]

    vxB = Vector{String}(undef,length(xyc))
    vyB = Vector{String}(undef,length(xyc))

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
