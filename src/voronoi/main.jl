using VoronoiDelaunay, Gadfly
include("vorfun.jl")
include("../geomlibs.jl")


xy = 100*rand(16,2);
xy0 = collect(Iterators.product(0:10:99,0:10:99))[:];
  xy = zeros(Float64,length(xy0),2);
  for (k,v) in enumerate(xy0) xy[k,:].=v; end

bnd  = [0 0; 100 0; 100 100; 0 100; 0 0]
tess, rc, Cv1, Cv = makeCell(xy,bnd)

x, y = getplotxy(Cv1)
x, y = getplotxy(Cv)
xt, yt = getplotxy(delaunayedges(tess))


set_default_plot_size(15cm, 15cm)
plot(x=x, y=y, Geom.path, Coord.cartesian(xmin=0, xmax=3, ymin=0, ymax=3))

plot(layer(x=x, y=y, Geom.path),
     layer(x=xt, y=yt, Geom.path,Theme(default_color="orange")),
     #layer(x=[Cv1[42]._generator_a._x,Cv1[42]._generator_b._x],
     #y=[Cv1[42]._generator_a._y,Cv1[42]._generator_b._y]),
     #layer(x=map(x->x[1],xyc),y=map(x->x[2],xyc)),
     layer(x=map(x->x[1],xyc)[bc],y=map(x->x[2],xyc)[bc]),
     Coord.cartesian(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5))

plot(layer(x=x, y=y, Geom.path),
  layer(x=xt, y=yt, Geom.path,Theme(default_color="orange")),
  layer(x=xy[:,1]/100 .+1, y=xy[:,2]/100 .+1),
  Coord.cartesian(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5))

plot(x=bnd[:,1], y=bnd[:,2], Geom.path)

bc = falses(length(xyc))
for i=1:length(xyc)
    fl.=false;
    ir = findall((xca.==xyc[i][1]) .& (yca.==xyc[i][2]))
    ic = findall((xcb.==xyc[i][1]) .& (ycb.==xyc[i][2]))
    fl[ir].=true;
    fl[ic].=true;
    #println(sum(fl))
    if any(bnd_flag[fl,1].!=0)
        ia = findall(bnd_flag[fl,1].!=0)
        println("$i _$(length(unique(bnd_flag[fl,1][ia])))")
        bc[i]= true
    end

end
