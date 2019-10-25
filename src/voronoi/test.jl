function test()

    xy = 100*rand(1024,2);
    bnd  = [0 0; 100 0; 100 100; 70 100; 70 70; 30 70; 30 100; 0 100; 0 0]

    tess, rc, Cv, exy = test_grid(xy,bnd);
    test_export(xy, rc, Cv, exy);
    println("all test is OK")
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
    return true
end

test()
