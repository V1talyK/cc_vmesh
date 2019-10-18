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
    Cv1, bnd_flag = makeCellvsBondary(Cv,bx,by);

    rc = makeRC(Cv1,x,y)

    return tess, rc, Cv1, Cv
end
function getPointToCell(Cv,xc,yc)
    xca = map(v->v._generator_a._x,Cv)
    xcb = map(v->v._generator_b._x,Cv)
    yca = map(v->v._generator_a._y,Cv)
    ycb = map(v->v._generator_b._y,Cv)

    xa = map(v->v._a._x,Cv)
    xb = map(v->v._b._x,Cv)
    ya = map(v->v._a._y,Cv)
    yb = map(v->v._b._y,Cv)

    xcf = vcat(xca,xcb)
    ycf = vcat(yca,ycb)
    xyc = unique(collect(zip(xcf,ycf)))


    fl = falses(length(xca))
    pnts = Vector(undef,length(xc))
    for i=1:length(xyc)
        fl.=false;
        ir = findall((xca.==xyc[i][1]) .& (yca.-yc[i].<1e-10))
        ic = findall((xcb.==xyc[i][2].<1e-10) .& (ycb.-yc[i].<1e-10))
        fl[ir].=true;
        fl[ic].=true;

        pnts[i] = [xa[fl],xb[fl]],[ya[fl],yb[fl]]
        if any(bnd_flag[fl].!=0)
            ia = findall(bnd_flag[fl].!=0)
            if length(unique(bnd_flag[fl][ia]))>1
                px1 = xa[fl][ia]
                ibo = maximum(unique(bnd_flag[fl][ia]))
                bxy[ibbo[ibo],:]
            else
                px1 = xa[fl][ia]
            end
        end

        if (length(ir)>0) .& (length(ic)>0)
            k+=1;
            rc[k,:] = [ir[1],ic[1]]
        end
    end

    return
end
function makeCellvsBondary(Cv,bx,by)
    Cv, ia, xy_ab, ab = remPointOutBondary(Cv,bx,by);
    ABC_l1 = lineEq(xy_ab[:,1:2],xy_ab[:,3:4])

    bxy = hcat(bx,by);
    ABC_l2 = lineEq(view(bxy,1:length(bx)-1,:),view(bxy,2:length(bx),:))

    bnd_flag = zeros(Int32,length(Cv))

    for j=1:length(bx)-1
        xy2 = Float64.([bx[j],by[j],bx[j+1],by[j+1]])
        flag_int = SegIntersect(xy_ab,xy2)
        ib = findall(flag_int)
        bnd_flag[findall(ia)[ib]] .= j;
        for (k,v) in enumerate(flag_int)
            if v
                xi = lineIntersect(ABC_l1[k,:]',ABC_l2[j,:]')
                if ab[k]
                    xy_ab[k,1:2]=xi
                else
                    xy_ab[k,3:4]=xi
                end
            end
        end
    end

    Cv1=copy(Cv)
    k1=0;
    for (k,v) in enumerate(ia)
        if v
          k1+=1;
          Cv1[k]=VoronoiDelaunay.VoronoiEdge(
                    Point2D(xy_ab[k1,1], xy_ab[k1,2]),
                    Point2D(xy_ab[k1,3], xy_ab[k1,4]),
                    Point2D(Cv[k]._generator_a._x, Cv[k]._generator_a._y),
                    Point2D(Cv[k]._generator_b._x, Cv[k]._generator_b._y))
        end
    end
    return Cv1, bnd_flag
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
function remPointOutBondary(Cv,bx,by)
    xa = map(v->v._a._x,Cv)
    xb = map(v->v._b._x,Cv)
    ya = map(v->v._a._y,Cv)
    yb = map(v->v._b._y,Cv)

    pin1 = inPolygon(xa, ya, bx, by)
    pin2 = inPolygon(xb, yb, bx, by)
    ic = .!(.!pin1 .& .!pin2);
    ia = .!(pin1 .& pin2)[ic]

    ab = falses(length(pin1));
    ab[pin1.==false] .=true;
    ab[pin2.==false] .=false;

    xy = hcat(xa,ya,xb,yb)
   return Cv[ic], ia, xy[ic,:][ia,:], ab[ic][ia]
end
