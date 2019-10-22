function make_vfile(OUT,rc,Cv)
    xyc_ab, xy_ab = getFromCV(Cv1);
    ic_ab, xyc = get_index(xyc_ab); #индексы центров и координаты
    ie_ab, xye = get_index(xy_ab);  #индексы границ и координаты

    i_bnd = findall((x->all(x==(0.,0.))).(xyc))[1]; #Индекс отвечающий за границу
    i_ed = (x->x[1]).(findall(ic_ab.==i_bnd)); #Индексы связей, которые являются границей;

    bnd_ind = zeros(Int32,length(i_ed),2)
    ebo = eage_order(xye,ie_ab,i_bnd,ic_ab); #порядок индексов границы иджей из xye
    indexin(ie_ab,ebo)


    bnd_ind[1,1] = setdiff(ic_ab[ia[1],:],i_bnd)[1];

    ib = findall(any(bnd_ind[1,1].==ic_ab,dims=2)[:] .& any(ind_1.==ie_ab,dims=2)[:])
    ind_2 = setdiff(ie_ab[ib[1],:],ind_1)[1]

    ia = (x->x[1]).(findall(ie_ab.==ind_2));
    bnd_ind[1,2] = setdiff(unique(ic_ab[ia,:]),[i_bnd,bnd_ind[1,1]])[1];
    bnd_ind[2,1]=bnd_ind[1,2]

    for i=2:6
        println(i)
        global ind_2
        ind_1 = copy(ind_2)
        ib = findall(any(bnd_ind[i,1].==ic_ab,dims=2)[:] .& any(ind_1.==ie_ab,dims=2)[:])

        ind_2 = setdiff(ie_ab[ib[1],:],ind_1)[1]

        ia = (x->x[1]).(findall(ie_ab.==ind_2));
        bnd_ind[i,2] = setdiff(unique(ic_ab[ia,:]),[i_bnd,bnd_ind[i,1]])[1];
        bnd_ind[i+1,1]=bnd_ind[i,2]
    end


    ia = indexin(xyc,xy[bnd_ind]);
    ia = reshape(ia,length(Cv1),2)
    ia = any(ia.!=nothing,dims=2)[:]

    num_bo = zeros(Int64,length(bnd_ind));
    num_bo1 = Vector(undef,0);

    k1=0
    for (k,v) in enumerate(ia)
        if v
            if all(xyc_ab[k,2:2:4].==0.)
                global k1+=1
                push!(num_bo1,[k, xy_ab[k,:]])
            end
        end
    end

    id1 = findall((x->any(x[2][1:2].==1) & any(x[2][3:4].==1)).(num_bo1))[1]
    next_p = [view(num_bo1[id1][2],1:2:3),view(num_bo1[id1][2],2:2:4)].==[1.,1.]
    fl = true
    while fl

    end

    make_mesh(xy,bnd_ind,OUT[1])
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

function get_index(xy_ab)
    xyab = collect(zip(view(xy_ab,:,1:2)[:],view(xy_ab,:,3:4)[:]))
    xy = unique(xyab); #xy =xy[.!(x->all(x==(0.,0.))).(xy)]
    ia = indexin(xyab,xy)
    return reshape(ia,size(xy_ab,1),2), xy
end

function make_mesh(xy,OUT)
    n = length(xy)
    id = collect(1:length(xy))
    #xy =  convert(Array{Tuple{Float64,Float64},2},xy)
    bnd = zeros(Int32,size(xy));
    # id50 = reshape(id,50,50)
    # bnd[1,1:end-1] .= id50[1,2:end];
    # bnd[1:end-1,end] .= id50[2:end,end];
    # bnd[end,2:end] .= id50[end,1:end-1];
    # bnd[2:end,1] .= id50[1:end-1,1];

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
