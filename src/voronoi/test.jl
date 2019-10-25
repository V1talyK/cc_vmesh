function test(v=1,b=1; remove::Bool=true)
    if v==1
        xy = 100*rand(1024,2);
    else
        #Прямоугольная сетка
        xy0 = collect(Iterators.product(0:10:99,0:10:99))[:];
          xy = zeros(Float64,length(xy0),2);
          for (k,v) in enumerate(xy0) xy[k,:].=v; end
    end
    if b==1
        bnd  = [0 0; 100 0; 100 100; 0 100; 0 0]
    elseif b==2
        bnd  = [0 0; 100 0; 100 70; 70 100; 0 100; 0 0]
    else
        bnd  = [0 0; 100 0; 100 100; 70 100; 70 70; 30 70; 30 100; 0 100; 0 0]
    end

    tess, rc, Cv, exy = test_grid(xy,bnd);
    OUT = test_export(xy, rc, Cv, exy);
    println("all test is OK")
    if remove
        map(x->rm(x),OUT)
    end

    return true
end

function test_grid(xy,bnd)
    @time tess, rc, Cv, Cv0, exy = makeCell(xy,bnd)
    println("test of makeVoronoi is OK")
    return tess, rc, Cv, exy
end

function test_export(xy,rc,Cv,exy)
    r = dirname(Base.source_path());
    OUT = [joinpath(r,"4_test_mesh.tsv"),joinpath(r,"5_test_geom.tsv"),joinpath(r,"6_test_wellCon.tsv")]
    wxy = collect(Iterators.partition(xy',2))
    make_vfile(OUT,wxy,rc,Cv, exy)
    println("test of export is OK")
    return OUT
end

test()
