function makeCell(xy)
    exy = extrema(xy,dims=1)
    x = 1 .+(xy[:,1].-exy[1][1])./exy[1][2];
    y = 1 .+(xy[:,2].-exy[2][1])./exy[2][2];
    x = clamp.(x,1+eps(Float64),2-2*eps(Float64))
    y = clamp.(y,1+eps(Float64),2-2*eps(Float64))

    tess = DelaunayTessellation(length(x))
    a = Point2D[Point(i[1], i[2]) for i in zip(x,y)]
    @time for i=1:length(a)
        println(i)
        push!(tess, a[i])
    end

    chan = voronoiedges(tess)

    Cv = [i for i in chan]
    # x1 = map(v->v._a._x,Cv)
    # x2 = map(v->v._b._x,Cv)
    # y1 = map(v->v._a._y,Cv)
    # y2 = map(v->v._b._y,Cv)
    rc = makeRC(Cv,x,y)

    return tess, rc
end

function makeRC(Cv,xc,yc)
    xg1 = map(v->v._generator_a._x,Cv)
    xg2 = map(v->v._generator_b._x,Cv)
    yg1 = map(v->v._generator_a._y,Cv)
    yg2 = map(v->v._generator_b._y,Cv)

    rc = zeros(Int64,length(xg1),2)#Array{Int64,2}(undef,length(xg1),2)
    k=0
    for i=1:length(xg1)
        ir = findall((xg1[i].==xc) .& (yg1[i].==yc))
        ic = findall((xg2[i].==xc) .& (yg2[i].==yc))
        if (length(ir)>0) .& (length(ic)>0)
            k+=1;
            rc[k,:] = [ir[1],ic[1]]
        end
    end
    return rc[1:k,:]
end
