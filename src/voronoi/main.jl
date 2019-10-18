using VoronoiDelaunay, Gadfly
include("vorfun.jl")
include("../geomlibs.jl")

plot(x=xg1, y=yg1, Coord.cartesian(xmin=1, xmax=2, ymin=1, ymax=2))
plot(x=map(x->x._x,a), y=map(x->x._y,a), Coord.cartesian(xmin=1, xmax=2, ymin=1, ymax=2))


VC = Vector(undef,100)
for i=1:100
    ia = findall((xc[i].==xg1) .& (yc[i].==yg1))
    VC[i] = (x=[x1[ia] x2[ia]], y = [y1[ia] y2[ia]])
end

xy = 100*rand(16,2);
xy0 = collect(Iterators.product(0:10:99,0:10:99))[:];
  xy = zeros(Float64,length(xy0),2);
  for (k,v) in enumerate(xy0) xy[k,:].=v; end

bnd  = [0 0; 100 0; 100 100; 0 100; 0 0]
tess, rc, Cv1, Cv = makeCell(xy,bnd)

x, y = getplotxy(Cv)
xt, yt = getplotxy(delaunayedges(tess))


set_default_plot_size(15cm, 15cm)
plot(x=x, y=y, Geom.path, Coord.cartesian(xmin=0, xmax=3, ymin=0, ymax=3))

plot(layer(x=x, y=y, Geom.path),
     layer(x=xt, y=yt, Geom.path,Theme(default_color="orange")),
     #layer(x=[Cv1[42]._generator_a._x,Cv1[42]._generator_b._x],
     #y=[Cv1[42]._generator_a._y,Cv1[42]._generator_b._y]),
     layer(x=map(x->x[1],xycp),y=map(x->x[2],xycp)),
     Coord.cartesian(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5))

plot(layer(x=x, y=y, Geom.path),
  layer(x=xt, y=yt, Geom.path,Theme(default_color="orange")),
  layer(x=xy[:,1]/100 .+1, y=xy[:,2]/100 .+1),
  Coord.cartesian(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5))

plot(x=bnd[:,1], y=bnd[:,2], Geom.path)

bc = falses(length(xcf))
for i=1:length(xcf)
    fl.=false;
    ir = findall((xca.==xcf[i]) .& (yca.==ycf[i]))
    ic = findall((xcb.==xcf[i]) .& (ycb.==ycf[i]))
    fl[ir].=true;
    fl[ic].=true;
    println(sum(fl))
    if any(bnd_flag[fl].!=0)
        ia = findall(bnd_flag[fl].!=0)
        println("$k _$(length(unique(bnd_flag[fl][ia])))")
        bc[i]= true
    end

end
