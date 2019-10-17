function makeCell(xy, bnd)

    fl = inPolygon(view(xy,:,1), view(xy,:,2), view(bnd,:,1), view(bnd,:,2))
    xy = xy[fl,:];

    exy = extrema(bnd,dims=1)
    x = 1 .+(xy[:,1].-exy[1][1])./exy[1][2];
    y = 1 .+(xy[:,2].-exy[2][1])./exy[2][2];
    x = clamp.(x,1+eps(Float64),2-2*eps(Float64))
    y = clamp.(y,1+eps(Float64),2-2*eps(Float64))

    bx = 1 .+(bnd[:,1].-exy[1][1])./exy[1][2];
    by = 1 .+(bnd[:,2].-exy[2][1])./exy[2][2];

    tess = DelaunayTessellation(length(x))
    a = Point2D[Point(i[1], i[2]) for i in zip(x,y)]
    @time for i=1:length(a)
        println(i)
        push!(tess, a[i])
    end

    chan = voronoiedges(tess)

    Cv = [i for i in chan]
    x1 = map(v->v._a._x,Cv)
    x2 = map(v->v._b._x,Cv)
    y1 = map(v->v._a._y,Cv)
    y2 = map(v->v._b._y,Cv)

    xy1 = hcat(x1,y1,x2,y2)
    bxy = hcat(bx,by);
    ABC_l2 = lineEq(view(bxy,1:length(bx)-1,:),view(bxy,2:length(bx),:))

    pin1 = inPolygon(x1, y1, bx, by)
    pin2 = inPolygon(x2, y2, bx, by)
    ia = findall(.!(pin1 .& pin2))
    ic = .!(.!pin1 .& .!pin2);

    for j=1:length(bx)-1
        xy2 = Float64.([bx[j],by[j],bx[j+1],by[j+1]])
        flag_int = SegIntersect(xy1[ia,:],xy2)
        ib = findall(flag_int)

        ABC_l1 = lineEq(xy1[ia[ib],1:2],xy1[ia[ib],3:4])
        k=0
        for i=1:length(flag_int)
            if flag_int[i]
                k+=1;
                xi = lineIntersect(ABC_l1[k,:]',ABC_l2[j,:]')
                fl = inPolygon(xy1[ia[i],1:2:3], xy1[ia[i],2:2:4], bx, by)
                if !fl[1]
                    x1[ia[i]]=xi[1]
                    y1[ia[i]]=xi[2]
                elseif !fl[2]
                    x2[ia[i]]=xi[1]
                    y2[ia[i]]=xi[2]
                end
            end
        end
    end

    Cv1=Vector(undef, length(Cv))
    for (k,v) in enumerate(Cv)
        Cv1[k]=VoronoiDelaunay.VoronoiEdge(Point2D(x1[k], y1[k]),
                                      Point2D(x2[k], y2[k]),
                                      Point2D(v._generator_a._x, v._generator_a._y),
                                      Point2D(v._generator_b._x, v._generator_b._y))

        # v._a._x = x1[k]
        # v._b._x = x2[k]
        # v._a._y = y1[k]
        # v._b._y = y2[k]
    end

    rc = makeRC(Cv,x,y)

    return tess, rc, Cv1[ic], Cv
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
