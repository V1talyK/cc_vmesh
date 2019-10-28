using VoronoiDelaunay, Gadfly, Dates
include("vorfun.jl")
include("make_vfile.jl")
include("../geomlibs.jl")
include("../libs.jl")
include("test.jl")

r = joinpath(dirname(dirname(dirname(Base.source_path()))),"data")
OUT = [joinpath(r,"4_mesh.tsv"),joinpath(r,"5_geom.tsv"),joinpath(r,"6_wellCon.tsv")]

#xy = 100*rand(1024,2);

@time tess, rc, Cv, Cv0, exy = makeCell(wxy,bnd)
make_vfile(OUT,wi,wxy,rc,Cv, exy)


x, y = getplotxy(Cv)
x, y = getplotxy(Cv0)
xt, yt = getplotxy(delaunayedges(tess))


set_default_plot_size(15cm, 15cm)
plot(x=x, y=y, Geom.path, Coord.cartesian(xmin=1., xmax=2., ymin=1., ymax=2.))


plot(layer(x=x, y=y, Geom.path),
     layer(x=xt, y=yt, Geom.path,Theme(default_color="orange")),
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

plot(layer(x=(x->x[1]).(xyc)[pbo],y=(x->x[2]).(xyc)[pbo], Geom.path),
     layer(x=x, y=y, Geom.path,Theme(default_color="orange")))

plot(layer(x=(x->x[1]).(xyc)[pbo],y=(x->x[2]).(xyc)[pbo], Geom.path),
     layer(x=(x->x[1]).(wxy),y=(x->x[2]).(wxy), Geom.point),
     layer(x=mk_decmprs(x,exy[1]), y=mk_decmprs(y,exy[2]), Geom.path,Theme(default_color="orange")))

plot(layer(x=wxy[:,1],y=wxy[:,2]),
          layer(x=bnd[:,1], y=bnd[:,2], Geom.path,Theme(default_color="orange")))
