function makeVecOfDat(obj)
    ds= floor(Date(obj["date_start"]), Dates.Month);
    de= floor(Date(obj["date_end"]), Dates.Month);
    vd=ds:Dates.Month(1):de;
    return vd
end
function tmap(fun, args...)
  RSL = Vector(undef, length(args[1]))
  @inbounds Threads.@threads for i7 = 1:length(args[1])
     RSL[i7] = fun(map(x->x[i7],args)...)
  end
 return RSL
end
function date_iso(diso::Union{String,SubString{String}})
    #Ускоренное преобразование панорамных дат формата YYYY-MM-DD
    #println(diso)
    try
        year = parse(Int, diso[1:4])
        month = parse(Int, diso[6:7])
        day = length(diso)>=10 ? parse(Int, diso[9:10]) : 1;
        return Date(year, month, day)
    catch
        println(diso)
    end
end

function writeToPipe(io, args...)
    k = 0

    newline = '\n'
    delim = '\t'
    @inbounds while k < length(args[1])
        k += 1
        ncols = length(args)
        @inbounds for i = 1:ncols
            v = getStr(args[i][k])
            write(io, "$v")
            write(io, ifelse(i == ncols, newline, delim))
        end

    end
    close(io)
    return true
end

function writeToPipe(io, vd::StepRange, args::Array{Float64,2})
    k = 0

    newline = '\n'
    delim = '\t'
    @inbounds while k < length(vd)
        k += 1
        #ncols = length(args)
        #@inbounds for i = 1:ncols
            v = getStr(vd[k])
            write(io, "$v")
            write(io, delim)
            #v = getStr(args[k,:][:])
            v=string(args[k,:])
            v= replace(v, "[" => "\"")
            v = replace(v, "]" => "\"")
            write(io, v)
            write(io, newline)
        #end

    end
    close(io)
    return true
end

@inline getStr(v::Number) = ifelse(isnan(v), "\\N", string(v))
#@inline getVal(v::Array{T,1}) where T <: Number = (length(v) > 0) ? "["*join(map(x -> string(x), v),", ")*"]" : "\\N"
#@inline getVal(v::Array{T,1}) where T <: AbstractString = (length(v) > 0) ? "["*join(map(x -> "'$x'", v),", ")*"]" : "\\N"
@inline getStr(v::Tuple) = string(v)
@inline getStr(v::Date) = "$v"
@inline getStr(v::String) = "$v"
@inline getStr(v::Array{Float64,1}) = "$v"
@inline getStr(v::Array{<:Signed,1}) = "$v"
#@inline getVal(v::Nothing) = "\\N"
#@inline getVal(v::Missing) = "\\N"
#@inline getVal(v::Array{Missing,1}) = "\\N"

function readFromPipe(tmp,d_ind=2,f_ind=3)

    tmp = split(tmp);

    ncols = length(tmp)
    v = zeros(Float64,ncols-f_ind+1)
    wi = parse(Int32,tmp[1])
    d = date_iso.(tmp[d_ind])

    k = 0
    @inbounds for i = f_ind:ncols
        k+=1;
        z = tryparse(Float64,tmp[i]);
        v[k] = z!=nothing ? z : NaN
    end
    return wi, d, v
end

function readFromPipe2(tmp,args)
    k = 0
    tmp = split(tmp);
    ncols = length(tmp)
    v = zeros(Float64,ncols-1)
    wi = parse(Int32,tmp[1])
    @inbounds for i = 2:ncols
        z = tryparse(Float64,tmp[i]);
        v[i-1] = z!=nothing ? z : NaN
    end
    return wi, v
end

function readFromPipe3(tmp,args)
    k = 0
    tmp = split(tmp,"\t");
    ncols = length(tmp)
    v = Vector{Array{Real,1}}(undef,ncols-1)
    wi = parse(Int32,tmp[1])
    @inbounds for i = 2:ncols
        pt = Meta.parse(tmp[i]; raise=false);
        if isa(pt,Number)
            z = tryparse(Float64,tmp[i]);
            v[i-1] = z!=nothing ? [z] : [NaN]
        elseif isa(pt,String)
            v[i-1] = Meta.parse.(split(pt,","))
        else
            v[i-1] = [NaN];
        end
    end
    return wi, v
end

@inline getVal(x,t) = ifelse(t, "\\N", string(v))
#@inline getVal(v::Array{T,1}) where T <: Number = (length(v) > 0) ? "["*join(map(x -> string(x), v),", ")*"]" : "\\N"
#@inline getVal(v::Array{T,1}) where T <: AbstractString = (length(v) > 0) ? "["*join(map(x -> "'$x'", v),", ")*"]" : "\\N"
@inline getVal(v::Tuple) = string(v)
@inline getVal(v::Date) = "$v"
@inline getVal(v::String) = "$v"


function writeTSV(pipe::String, args...)

    io = open(pipe, "w")
    writeToPipe(io, args...)

end

function getInputCmd(s3ng::String)
    a = split(s3ng);
    i_fun = findfirst(a.=="--fun")
    fun = a[i_fun+1];
    fl_in = find(map(x->contains(x,"--in"),a));
    IN = getInOut(a[fl_in], a[fl_in+1],"in")

    fl_out = find(map(x->contains(x,"--out"),a));
    OUT = getInOut(a[fl_out], a[fl_out+1],"out")
    return fun, IN, OUT
end

function getInOut(k::Array{SubString{String},1},v::Array{SubString{String},1},s3ng)
    fl = map(x->contains(x,s3ng),k);
    ord = map(x->x[search(x,s3ng)[end]+1:end],k[fl]);
    ib = find(ord.!="");
    ic = find(ord.=="");
    ia = parse.(Int,ord[ib]);
    return vcat(v[fl][ic],v[fl][ia])
end

function CSVread(ioR,to=[Int32,Date,Float64])

    Wi=zeros(Int32,0);
    D=Vector{Date}(undef,0);
    V = Vector{Array{Float64,1}}(undef,0);
    d_ind = findall(to.==Date);
    f_ind = findfirst(to.==Float64)
    while !eof(ioR)
        tmp = Base.readline(ioR)
        wi, d, v = readFromPipe(tmp,d_ind,f_ind)
        push!(Wi,wi)
        push!(D,d[1])
        push!(V,v)
    end
    return Wi, D, V
end

function CSVread2(ioR)
    k = 0
    Wi=zeros(Int32,0);
    V = Vector{Array{Float64,1}}(undef,0);
    while !eof(ioR)
        k+=1;
        tmp = Base.readline(ioR)
        wi, v = readFromPipe2(tmp,[Int32,Date,Float64])
        push!(Wi,wi)
        push!(V,v)
    end
    return Wi, V
end

function CSVreadVec(ioR)
    k = 0
    Wi=zeros(Int32,0);
    V = Vector(undef,0);
    while !eof(ioR)
        k+=1;
        tmp = Base.readline(ioR)
        wi, v = readFromPipe3(tmp,[Int32,Date,Float64])
        push!(Wi,wi)
        push!(V,v)
    end
    return Wi, V
end

function CSVread3(ioR,to=[Int32,Date,Float64])

    Wi=zeros(Int32,0);
    D=Vector{Array{Date,1}}(undef,0);
    V = Vector{Array{Float64,1}}(undef,0);
    d_ind = findall(to.==Date);
    f_ind = findfirst(to.==Float64)
    while !eof(ioR)
        tmp = Base.readline(ioR)
        wi, d, v = readFromPipe(tmp,d_ind,f_ind)
        push!(Wi,wi)
        push!(D,d)
        push!(V,v)
    end
    return Wi, D, V
end
function dinFromTlb(tlb,clmTps::Array{DataType,1})
    lTlb = length(tlb[1]);

    wi = tlb[1];
    d = tlb[2];
    DV = tlb[3];
    DV = hcat(DV...)';

    ds= floor(minimum(d), Dates.Month);
    de= floor(maximum(d), Dates.Month);
    vd=ds:Dates.Month(1):de;

    uwi = unique(wi);

   nw=length(uwi)
   nt=length(vd);

  r = indexin(d,vd);
  c = indexin(wi,uwi);
  #pw, p, pdl, psl, pb, uf, qw, qo = map(x->full(sparse(r,c,x,nt,nw)),DV);

 return uwi, vd, map(x->Array(sparse(r,c,DV[:,x],nt,nw)),1:size(DV,2))
end
function statFromTlb(tlb,clmTps::Array{DataType,1})
    wi = tlb[1];
    DV=hcat(tlb[2]...)'
  return wi, DV
end
function vecFromTlb(tlb,clmTps::Array{DataType,1})
    wi = tlb[1];
    DV = Vector(undef,length(clmTps)-1)
    for i=1:length(clmTps)-1
        DV[i]=map(x->x[i],tlb[2])
    end
  return wi, DV
end

function getFMod(res)
    return getFMod(res["fmod_id"],res["pvt"])
end

function getFMod(fmod_id,pvt)
    fmod_id=Int16(fmod_id);
    fmod_prm=getFLUIDPrm(pvt);
    return fmod_id, fmod_prm
end

function getFLUIDPrm(fmod_prm)
    #fl=Vector{SLX.FL}(undef,3)
    Ro = Base.get(fmod_prm,"Density",Dict("oil"=>900,"water"=>1000,"gaz"=>0.656))  #Плотность
    Mu = Base.get(fmod_prm,"Viscosity",Dict("oil"=>1,"water"=>1,"gaz"=>1.1e-5))   #Вязкоасть
    Comp = Base.get(fmod_prm,"Compressibility",Dict("oil"=>1e-5,"water"=>1e-6,"gaz"=>1e-3)) #Сжимаемость

    fl1 = FL(Ro["water"],Ro["oil"],Ro["gaz"],NaN)
    fl2 = FL(Mu["water"],Mu["oil"],Mu["gaz"],NaN)
    fl3 = FL(Comp["water"],Comp["oil"],Comp["gaz"],NaN)

    return FLPRM(fl1,fl2,fl3)
end

function getDataFromAdp(obj)
  APC=Base.get(obj,"AdpPrm",obj)
  PRM = Dict();
  PRM["maxIterP"] = Base.get(APC, "maxIterP",50);
  PRM["maxIterG"] = Base.get(APC, "maxIterG",100);
  PRM["dpw"] = Base.get(APC, "dPw",15)/100;
  PRM["dpw2"] = Base.get(APC, "dPw2",20)/100;
  PRM["dpw3"] = Base.get(APC, "dPw3",25)/100; #Аварийный коэффициент, когда упёрлось в ограничение - неверный замер
  PRM["dqp"] = Base.get(APC, "dQp",10)/100;
  PRM["dqi"] = Base.get(APC, "dQi",10)/100;
  PRM["kqi"] = Base.get(APC, "K_Qi",1); #Допустимый коэффициент потерь закачки
  PRM["maxG"] = 10e12; #Максимальная гидропроводность
  PRM["skinFlag"] = Base.get(APC, "useSkin",false); #Адаптировать скинфактор
  PRM["fl_get_slx"] = Base.get(APC, "fl_get_slx",false); #Использовать предварительный расчёт гидропроводности
  PRM["piezFlag"] = Base.get(APC, "usePiezo",false); #Использовать предварительный расчёт гидропроводности
  #PRM=[APC["max_k"],APC["dPw"]/100,APC["dPw2"]/100, 0.5,APC["dQp"]/100,APC["dQi"]/100]
  PRM["snt"] = Base.get(APC, "set_num_threads",false); #Изменить дефолтное кол-во потоков
  PRM["num_threads"] = Base.get(APC, "num_threads",1); #Кол-во потоков при последовательном расчёте

  PRM["useBottomQ"] = Base.get(APC, "useBottomQ",false); #Приток через подошву
  PRM["partBottomQ"] = Base.get(APC, "partBottomQ",0); #Доля притока через подошву

  return PRM
end

function checkFM(id,prm,pcol...)
    p=Vector(undef,length(pcol)-1);
    flag=pcol[end];
    if flag==1 #конвертим в нужную функцию
        if id==Int16(1)
            p=pcol[1:end-1];
        elseif id==Int16(2)
            rog=prm.ro.g
            Pat=1;
            for i=1:length(pcol)-1
                p[i]=rog/Pat*pcol[i].^2/2
            end
        elseif id==Int16(4)
            # WNF=1;
            # E=GNF.*mug./muo;
            # H=Hz*Pat*E;
            p=pcol[1:end-1];
        elseif id==Int16(5)
            p=pcol[1:end-1];
        else
        end
    else #конвертим обратно в давление
        if id==Int16(1)
            p=pcol[1:end-1];
        elseif id==Int16(2)
            rog=prm.ro.g
            Pat=1;
            for i=1:length(pcol)-1
                p[i]=pcol[i];
                p[i][pcol[i].<0]=0;
                p[i]=sqrt(2*pcol[i]*Pat/rog)
            end
        elseif id==Int16(4)
            # WNF=1;
            # E=GNF.*mug./muo;
            # H=Hz*Pat*E;
            p=pcol[1:end-1];
        elseif id==Int16(5)
            p=pcol[1:end-1];
        else
        end
    end
    return p
end
function doPermut(grd)
    Ar = rand(length(grd.B))
    A=sparse(view(grd.rc,:,1),view(grd.rc,:,2),Ar);
    p = symrcm(A)
    #AP = rcmpermute(A)

    r = Int64.(indexin(grd.rc[:,1],p))
    c = Int64.(indexin(grd.rc[:,2],p))
    grd.Won[:,1] = Int32.(indexin(grd.Won[:,1],p))
    grd.XY = grd.XY[p,:]
    grd.H = grd.H[p]
    grd.Sp = grd.Sp[p]
    grd.rc = hcat(r,c);
    grd.sb = grd.sb[p]
    grd.p = p
   return grd
end

function nanmean(x)
    iv = .!isnan.(x);
    if any(iv)
        x=mean(view(x,iv))
        return x
    else
        return NaN
    end
end
# sum(abs.(A[p,p]-AP))
# sum(abs.(A1-AP))
# sum(abs.(sparse(r,c,grd.B)-B1[p,p]))
