using VoronoiDelaunay, Gadfly
include("vorfun.jl")

plot(x=xg1, y=yg1, Coord.cartesian(xmin=1, xmax=2, ymin=1, ymax=2))
plot(x=map(x->x._x,a), y=map(x->x._y,a), Coord.cartesian(xmin=1, xmax=2, ymin=1, ymax=2))


VC = Vector(undef,100)
for i=1:100
    ia = findall((xc[i].==xg1) .& (yc[i].==yg1))
    VC[i] = (x=[x1[ia] x2[ia]], y = [y1[ia] y2[ia]])
end

xy = 100*rand(10000,2);
xy0 = collect(Iterators.product(0:10:99,0:10:99))[:];
  xy = zeros(Float64,length(xy0),2);
  for (k,v) in enumerate(xy0) xy[k,:].=v; end

tess, rc = makeCell(xy)

x, y = getplotxy(voronoiedges(tess))
set_default_plot_size(15cm, 15cm)
plot(x=x, y=y, Geom.path, Coord.cartesian(xmin=1, xmax=2, ymin=1, ymax=2))
