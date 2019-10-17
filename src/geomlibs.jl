function polyaria(XY)
  S=abs(sum((XY[1:end-1,1]+XY[2:end,1]).*(XY[1:end-1,2]-XY[2:end,2])))/2;
  return S
end

function perim(XY)
  P=sum(Distances.colwise(Distances.Euclidean(),XY[1:end-1,:]',XY[2:end,:]'))
  return P
end

function pointOnLine(xy1,xy2,abc,xy)
  #Принадлежит ли точка отрезку записанного прямой Ax+By+C=0;
  # xy1, xy2 - начало и конец отрезка; xy - точка
  i=1;
  ab=1 ./(abc[1]^2+abc[2]^2);
  xyp=zeros(Float32,size(xy,1),2);
  xyp[:,1]=(abc[2].*(abc[2]*xy[:,1]-abc[1]*xy[:,2])-abc[1]*abc[3])*ab;
  xyp[:,2]=(abc[1].*(-abc[2]*xy[:,1]+abc[1]*xy[:,2])-abc[2]*abc[3])*ab;

  pon=all(((xyp.>xy1') & (xyp.<xy2')) | ((xyp.<xy1') & (xyp.>xy2')),2);
  ll=Distances.pairwise(Distances.Euclidean(),xyp[pon[:].==false,:]',[xy1 xy2]);
  v,ind = findmin(ll,2);
  r,c =ind2sub(size(xyp[pon[:].==false,:]),ind[:]);
  xyp[pon[:].==false,:]=[xy1 xy2]'[c,:];

  return xyp, pon
end

function SegIntersect(xy1::Array{Float64,2},xy2::Array{Float64,2})
    #Пересекаются ли отрезки
    flag = falses(size(xy1));
    # ab1 = zeros(Float64,size(xy1,1));
    # ab2 = zeros(Float64,size(xy1,1));
    # a = zeros(Float64,size(xy1,1),2);
    # b = zeros(Float64,size(xy1,1),2);
    # c = zeros(Float64,size(xy1,1),2);

    A = view(xy1,:,1:2);
    B = view(xy1,:,3:4);
    C = view(xy2,:,1:2);
    D = view(xy2,:,3:4);

    #fl1 = kosProd!(ab1, ab2, a, b, c, A,B,C,D);
    fl1, ab1, ab2 = kosProd(A,B,C,D);

    flag[:,1]=ab1.==0;
    flag[:,2]=ab2.==0;

    fl2, ab1, ab2 = kosProd(C,D,A,B);

    flag[:,3]=ab1.==0;
    flag[:,4]=ab2.==0;

    flag2 = falses(size(xy1))
    flag2[:,1] = skalProd(C,A,B).<=0;
    flag2[:,2] = skalProd(D,A,B).<=0;
    flag2[:,3] = skalProd(A,C,D).<=0;
    flag2[:,4] = skalProd(B,C,D).<=0;
    return (fl1 & fl2) | any(flag & flag2, 2)[:]
end

function SegIntersect(xy1::Array{Float64,2},xy2::Array{Float64,1})
    #Пересекаются ли отрезки
    flag = falses(size(xy1));
    # ab1 = zeros(Float64,size(xy1,1));
    # ab2 = zeros(Float64,size(xy1,1));
    # a = zeros(Float64,size(xy1,1),2);
    # b = zeros(Float64,size(xy1,1),2);
    # c = zeros(Float64,size(xy1,1),2);

    A = view(xy1,:,1:2);
    B = view(xy1,:,3:4);
    C = view(xy2,1:2);
    D = view(xy2,3:4);

    #fl1 = kosProd!(ab1, ab2, a, b, c, A,B,C,D);
    fl1, ab1, ab2 = kosProd(A,B,C,D);

    flag[:,1]=ab1.==0;
    flag[:,2]=ab2.==0;

    fl2, ab1, ab2 = kosProd(C,D,A,B);

    flag[:,3]=ab1.==0;
    flag[:,4]=ab2.==0;

    flag2 = falses(size(xy1))
    flag2[:,1] = skalProd(C,A,B).<=0;
    flag2[:,2] = skalProd(D,A,B).<=0;
    flag2[:,3] = skalProd(A,C,D).<=0;
    flag2[:,4] = skalProd(B,C,D).<=0;
    return (fl1 & fl2) | any(flag & flag2, 2)[:]
end

function SegIntersectSL(A::Array{Float64,2},B::Array{Float64,2},xy2::Array{Float64,1},
                        ab1, ab2)
    #Пересекаются ли отрезки, кастом для линий тока
    flag = falses(size(A,1),4);
    flag2 = falses(size(A,1),4)
    C = view(xy2,1:2);
    D = view(xy2,3:4);

    fl1 = kosProd!(A,B,C,D,ab1, ab2);

    flag[:,1]=ab1.==0;
    flag[:,2]=ab2.==0;

    ia = findall(fl1);
    fl2=falses(length(fl1));
    fl2[ia] = kosProd!(C,D,view(A,ia,:),view(B,ia,:),view(ab1,ia), view(ab2,ia));

    flag[:,3]=ab1.==0;
    flag[:,4]=ab2.==0;

    ia = findall(flag[:,1]);
    flag2[ia,1] = skalProd(C,view(A,ia,:),view(B,ia,:)).<=0;
    ia = findall(flag[:,2]);
    flag2[ia,2] = skalProd(D,view(A,ia,:),view(B,ia,:)).<=0;
    ia = findall(flag[:,3]);
    flag2[ia,3] = skalProd(view(A,ia,:),C,D).<=0;
    ia = findall(flag[:,4]);
    flag2[ia,4] = skalProd(view(B,ia,:),C,D).<=0;
    return (fl1 .& fl2) .| any(flag .& flag2, dims = 2)[:]
end

function kosProd(A,B,C,D)
    a=B-A;
    b=C-A;
    c=D-A;

    ab1 = view(a,:,1).*view(c,:,2) - view(c,:,1).*view(a,:,2);
    ab2 = view(a,:,1).*view(b,:,2) - view(b,:,1).*view(a,:,2);
    fl1 = ab1 .* ab2.<0;
    return fl1, ab1, ab2
end
function kosProd(A,B,C::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                     D::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true})
    a=B-A;
    b=C'.-A;
    c=D'.-A;

    ab1 = view(a,:,1).*view(c,:,2) - view(c,:,1).*view(a,:,2);
    ab2 = view(a,:,1).*view(b,:,2) - view(b,:,1).*view(a,:,2);
    fl1 = ab1 .* ab2.<0;
    return fl1, ab1, ab2
end
function kosProd(A::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                 B::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                 C,D)
    a=B-A;
    b=C.-A';
    c=D.-A';

    ab1 = view(a,1).*view(c,:,2) - view(c,:,1).*view(a,2);
    ab2 = view(a,1).*view(b,:,2) - view(b,:,1).*view(a,2);
    fl1 = ab1 .* ab2.<0;
    return fl1, ab1, ab2
end
function kosProd!(A::Array{Float64,2},B::Array{Float64,2},
                  C::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                  D::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                  ab1,ab2)
    a=B-A;
    b=C'.-A;
    c=D'.-A;
    ab1[:] = view(a,:,1).*view(c,:,2) - view(c,:,1).*view(a,:,2);
    ab2[:] = view(a,:,1).*view(b,:,2) - view(b,:,1).*view(a,:,2);
    fl1 = ab1 .* ab2.<0;
    return fl1
end
function kosProd!(A::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                  B::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                  C::Union{SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Colon},false},
                    SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Base.Slice{Base.OneTo{Int64}}},false}},
                  D::Union{SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Colon},false},
                    SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Base.Slice{Base.OneTo{Int64}}},false}},
                  ab1,ab2)
    a=B-A;
    b=C.-A';
    c=D.-A';
    ab1[:] = view(a,1).*view(c,:,2) - view(c,:,1).*view(a,2);
    ab2[:] = view(a,1).*view(b,:,2) - view(b,:,1).*view(a,2);
    fl1 = ab1 .* ab2.<0;
    return fl1
end

function skalProd(A,B,C)
    a = B-A;
    b = C-A;
    return sum(a.*b,dims = 2)
end

function skalProd(A::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                  B,C)
    a = B.-A';
    b = C.-A';
    return sum(a.*b,dims = 2)
end

function skalProd(A,B::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                    C::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true})
    a = B'.-A;
    b = C'.-A;
    return sum(a.*b,dims = 2)
end

function lineEq(xy1,xy2)
# Рассчитывает коэффициенты уравнения прямой вида Ax+By+C=0
  ABC=zeros(Float64,size(xy1,1),3);
  ABC[:,1]=xy1[:,2]-xy2[:,2];
  ABC[:,2]=xy2[:,1]-xy1[:,1];
  ABC[:,3]=xy1[:,1].*xy2[:,2]-xy2[:,1].*xy1[:,2];
 return ABC
end

function lineIntersect(ABC1,ABC2)
    #Функция пересечения 2-х прямых записанных в виде Ax+By+C=0
    A1, B1, C1 = map(x-> view(ABC1,:,x), [1,2,3]);
    A2, B2, C2 = map(x-> view(ABC2,:,x), [1,2,3]);

    Z=(A1.*B2-A2.*B1);
    x=-(C1.*B2-C2.*B1)./Z;
    y=-(A1.*C2-A2.*C1)./Z;
    x[Z.==0].=NaN;
    y[Z.==0].=NaN;
    return hcat(x,y)
end

function pointOnSegment(xy1,xy2,xy)
  #Принадлежит ли точка отрезку
  # xy1, xy2 - начало и конец отрезка; xy - точка
  x, x1, x2 = map(z-> view(z,:,1), [xy, xy1, xy2]);
  y, y1, y2 = map(z-> view(z,:,2), [xy, xy1, xy2]);
  p1 = abs.((x-x1)./(x2-x1)-(y-y1)./(y2-y1)).<1e-5;
  p2 = ((x-x1).*(x-x2)+(y-y1).*(y-y2)).<=0
  return p1 .& p2
end

function cart2pol(x1, x0=zeros(Float64,size(x1)))
    #x1 - массив координат
    #x0 - массив координат начальных точек
    x=x1-x0;
    rho = sqrt.(sum(x.^2,dims = 2))
    phi = atan.(x[:,2],x[:,1])
    return rho[:], phi
end

function pol2cart(rho, phi)
    x=rho.*cos.(phi);
    y=rho.*sin.(phi);
    return hcat(x, y)
end

function sepPoint(xy, a, b)
    #Сортировка точек XY относительно прямой заданной 2-мя точками а и b
    #c=falses(size(xy,1));
    xy=xy.-a';
    b=b-a;
    c=(b[1]*xy[:,2]-b[2]*xy[:,1]).>0;
    return c
end

function distToLine(xy, a, b)
    #Расстояние от xy до прямой заданной 2-мя точками а и b
    A=(b[1]-a[1])*(xy[:,2].-a[2])-(b[2]-a[2])*(xy[:,1].-a[1]);
    B=sqrt((a[1]-b[1]).^2+(a[2]-b[2]).^2);
    h=abs.(A./B);
    return h
end

function distToSeg(xy, a, b)
    #Расстояние от xy до прямой заданной 2-мя точками а и b
    h=zeros(Float32, size(xy,1))
    f1 = funDistToSeg(xy, a, b)
    f2 = funDistToSeg(xy, b, a)
    fl = f1 .| f2
    h[.!fl] = distToLine(xy[.!fl,:], a, b);
    ha = sqrt.(sum((xy[fl,:].-a').^2,dims = 2))
    hb = sqrt.(sum((xy[fl,:].-b').^2,dims = 2))
    h[fl]=minimum(hcat(ha,hb),dims = 2)
    return h, .!fl, .!f1, .!f2
end

function funDistToSeg(xy, a, b)
    xy=xy.-a';
    b=b-a;
    fl = sum(xy.*b',dims = 2).<0
    return fl[:]
end

#@iprofile begin
function inPolygon(x, y, xp, yp)
    #Попадание точки x, y в полигон xp, yp
    if (xp[1]!=xp[end]) & (yp[1]!=yp[end])
        xp=vcat(xp,xp[1]);
        yp=vcat(yp,yp[1]);
    end
    fl1=falses(length(x))
    fl2=falses(length(x))
    c=falses(length(x))
    for i=2:length(xp)
       fl1 = ((yp[i].<=y) .& (y.<yp[i-1])) .| ((yp[i-1].<=y) .& (y.<yp[i]))
       fl2 = (x.> (xp[i-1] - xp[i]).* (y .- yp[i])./ (yp[i-1] - yp[i]) .+ xp[i])
       ia = findall(fl1 .& fl2)
       c[ia] = trues(length(ia)) - c[ia]
    end
   return convert(Array{Bool,1},c)
end

function inPolygon3(x, y, xp, yp)
    #Попадание точки x, y в полигон xp, yp
    if (xp[1]!=xp[end]) & (yp[1]!=yp[end])
        xp=vcat(xp,xp[1]);
        yp=vcat(yp,yp[1]);
    end
    fl1=falses(length(x))
    fl2=falses(length(x))
    c=falses(length(x))
    @inbounds for i=2:length(xp)
       fl1[:] = ((yp[i].<=y) & (y.<yp[i-1])) | ((yp[i-1].<=y) & (y.<yp[i]))
       ib = findall(fl1);
       fl2[:] = false;
       fl2[ib] = (view(x,ib).> (xp[i-1] - xp[i]).* (view(y,ib) - yp[i])./ (yp[i-1] - yp[i]) + xp[i])
       ia = findall(fl1 & fl2)
       c[ia] = trues(length(ia)) - c[ia]
    end
   return convert(Array{Bool,1},c)
end

#@iprofile begin
function convhull(XY, gi)
    #Оконтуривание точек выпуклой оболочкой
    XYgi = view(XY,gi,2);
    ci=gi[findmin(XYgi)[2]];
    P=copy(gi)
    H=ci;
    P=vcat(setdiff(P, H),H)

    fl=true
    k=0;
    while fl && k<length(gi)
      k+=1;
      right = 1
      Pr = view(P,right);
      XYPr = view(XY,Pr,:);
      XYH = XY[H[end],:];
      @inbounds for i = 1:length(P)
        # Pr = view(P,right);
        # XYPr = view(XY,Pr,:);
        if rotate(XYH,XYPr,XY[P[i],:])<0
          right = i
          Pr = view(P,right);
          XYPr = view(XY,Pr,:);
        end
      end

      if P[right]==H[1]
        fl==false
      else
        H=vcat(H,P[right])
        P=setdiff(P,P[right])
      end
    end
    return vcat(H,H[1])
end
function convhull_old(XY, gi)
    #Оконтуривание точек выпуклой оболочкой
    XYgi = view(XY,gi,2);
    ci=gi[findmin(XYgi)[2]];
    P=copy(gi)
    H=ci;
    P=vcat(setdiff(P, H),H)

    fl=true
    k=0;
    while fl && k<length(gi)
      k+=1;
      right = 1
      Pr = view(P,right);
      XYPr = view(XY,Pr,:);

      for i = 1:length(P)
        Pr = view(P,right);
        XYPr = view(XY,Pr,:);
        if rotate(XY[H[end],:],XYPr,XY[P[i],:])<0
          right = i
          # Pr = view(P,right);
          # XYPr = view(XY,Pr,:);
        end
      end

      if P[right]==H[1]
        fl==false
      else
        H=vcat(H,P[right])
        P=setdiff(P,P[right])
      end
    end
    return vcat(H,H[1])
end
function rotate(A,B,C)
  return (B[1]-A[1])*(C[2]-B[2])-(B[2]-A[2])*(C[1]-B[1])
end
#end
# type Properties
#
#     r::Array{Int64,1} #
#     c::Array{Int64,1}
#     b::Array{Float64,1}
#     areaEdge::Array{Float64,1}
#     distance::Array{Float64,1}
#
#     area::Array{Float64,1}
#     h::Array{Float64,1}
#
#     depth::Dict{String, Any}
#     data::Data
#
#     function Properties(data::Data)
#         new(zeros(Int64,0), zeros(Int64,0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), Dict(), data)
#     end
#
# end
#
# function idw(xy::Array{Float64,2},
#              z::Array{Float64,1},
#              xyi::Array{Float64,2},
#              nw::Int64,
#              nq::Int64,
#              delta::Float64;
#              method = 1
#             )
#
#     kdtree = KDTree(xy')
#
#     listNw, distNw = knn(kdtree, xyi', nw, true)
#     listNq, distNq = knn(kdtree, xy', nq, true)
#
#     wi = Array{Array{Float64,1},1}(length(listNw))
#     @inbounds for i = 1:length(listNw)
#         wi[i] = Array{Float64,1}(nw)
#         @inbounds for j = 1:nw
#             wi[i][j] = ( (distNw[i][end] - distNw[i][j]) / (distNw[i][end]*distNw[i][j] + delta) )^2
#         end
#     end
#
#     winq = Array{Array{Float64,1},1}(length(listNq))
#     @inbounds for i = 1:length(listNq)
#         winq[i] = Array{Float64,1}(nq - 1)
#         @inbounds for j = 2:nq
#             winq[i][j-1] = ( (distNq[i][end] - distNq[i][j]) / (distNq[i][end]*distNq[i][j] + delta) )^2
#         end
#     end
#
#     Ainv = getMatrixLS(winq, xy, listNq, method)
#
#     zi = zeros(length(listNw), size(z, 2))
#     for dim = 1:size(z, 2)
#         k = 0
#         @inbounds for i in listNw
#             k += 1
#             localSum = 0.
#             @inbounds for j = 1:nw
#                 neighb = listNq[i[j]][2:end]
#                 if method == 0
#                     b = [sum(winq[i[j]].*z[neighb, dim])]
#                 elseif method == 1
#                     b = [sum(winq[i[j]].*z[neighb, dim]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb ,1]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb, 2])
#                         ]
#                 else
#                     b = [sum(winq[i[j]].*z[neighb, dim]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb ,1]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb, 2]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb, 1].*xy[neighb, 2]),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb, 1].^2),
#                          sum(winq[i[j]].*z[neighb, dim].*xy[neighb, 2].^2)
#                         ]
#                 end
#                 coeff = Ainv[i[j]]*b
#                 if method == 0
#                     localSum += wi[k][j]*coeff[1]
#                 elseif method == 1
#                     localSum += wi[k][j]*(coeff[1] + coeff[2]*xyi[k,1] + coeff[3]*xyi[k,2])
#                 else
#                     localSum += wi[k][j]*(coeff[1] + coeff[2]*xyi[k,1] + coeff[3]*xyi[k,2] +
#                                 coeff[4]*xyi[k,1]*xyi[k,2] + coeff[5]*xyi[k,1]^2 + coeff[6]*xyi[k,2]^2)
#                 end
#             end
#
#            zi[k, dim] = localSum / sum(wi[k])
#
#         end
#     end
#     for dim = 1:size(z, 2)
#         zi[minimum(z[:, dim]) .> zi[:, dim], dim] = minimum(z[:, dim])
#         zi[maximum(z[:, dim]) .< zi[:, dim], dim] = maximum(z[:, dim])
#     end
#     return zi
# end
#
# function getMatrixLS(winq, xy, listNq, method=1)
#
#
#     k = 0
#     Ainv = Array{Array{Float64,2},1}(length(listNq))
#     x = xy[:,1]
#     y = xy[:,2]
#     for i in listNq
#         k += 1
#         i = i[2:end]
#         if method == 0
#             A = [sum(winq[k])]
#         elseif method == 1
#             A = [
#                  sum(winq[k]) sum(winq[k].*x[i]) sum(winq[k].*y[i]);
#                  sum(winq[k].*x[i]) sum(winq[k].*x[i].^2) sum(winq[k].*x[i].*y[i]);
#                  sum(winq[k].*y[i]) sum(winq[k].*x[i].*y[i]) sum(winq[k].*y[i].^2)
#                 ]
#         else
#             A = [
#                  sum(winq[k]) sum(winq[k].*x[i]) sum(winq[k].*y[i]) sum(winq[k].*x[i].*y[i]) sum(winq[k].*x[i].^2) sum(winq[k].*y[i].^2);
#                  sum(winq[k].*x[i]) sum(winq[k].*x[i].^2) sum(winq[k].*x[i].*y[i]) sum(winq[k].*y[i].*x[i].^2) sum(winq[k].*x[i].^3) sum(winq[k].*x[i].*y[i].^2);
#                  sum(winq[k].*y[i]) sum(winq[k].*x[i].*y[i]) sum(winq[k].*y[i].^2) sum(winq[k].*x[i].*y[i].^2) sum(winq[k].*y[i].*x[i].^2) sum(winq[k].*y[i].^3);
#                  sum(winq[k].*x[i].*y[i]) sum(winq[k].*y[i].*x[i].^2) sum(winq[k].*x[i].*y[i].^2) sum(winq[k].*(x[i].*y[i]).^2) sum(winq[k].*y[i].*x[i].^3) sum(winq[k].*x[i].*y[i].^3);
#                  sum(winq[k].*x[i].^2) sum(winq[k].*x[i].^3) sum(winq[k].*y[i].*x[i].^2) sum(winq[k].*y[i].*x[i].^3) sum(winq[k].*x[i].^4) sum(winq[k].*(y[i].*x[i]).^2);
#                  sum(winq[k].*y[i].^2) sum(winq[k].*x[i].*y[i].^2) sum(winq[k].*y[i].^3) sum(winq[k].*x[i].*y[i].^3) sum(winq[k].*(x[i].*y[i]).^2) sum(winq[k].*y[i].^4)
#                 ]
#         end
#         try
#             Ainv[k] = inv(A)
#         catch
#             Ainv[k] = pinv(A)
#         end
#      end
#      return Ainv
# end
#
# function getNearestNeighbore(xy::Array{Float64,2},
#                              z::Array{Float64,1},
#                              xyi::Array{Float64,2}
#                             )
#
#     kdtree = KDTree(xy')
#     list, dist = knn(kdtree, xyi, 1, false)
#
#     zi = zeros(length(list), size(z, 2))
#     for dim = 1:size(z, 2)
#         for i = 1:length(list)
#             zi[i, dim] = z[list[i][1]]
#         end
#     end
#     return zi
# end
#
# function check(data::Data,
#                geom::Properties
#               )
#     if data.data["directional"] == nothing
#         bot = maximum(geom.depth["z2"])
#         top = minimum(geom.depth["z1"])
#
#         wi = unique(data.data["wi"])
#         directional = Array{Any, 1}(length(wi))
#         for i = 1:length(wi)
#             directional[i] = Dict()
#             directional[i]["well_id"] = wi[i]
#
#             wxy = data.data["wxy"][data.data["wi"] .== wi[i], :]
#             directional[i]["data"] = Vector(2)
#             if size(wxy, 1) == 1
#                 directional[i]["data"][1] = [wxy[1,1]; wxy[1,2]; -top; top]
#                 directional[i]["data"][2] = [wxy[1,1]; wxy[1,2]; -bot; bot]
#             else
#                 dr = top + ((wxy[2,1] - wxy[1,1])^2 + (wxy[2,2] - wxy[1,2])^2 + (bot - top)^2)^0.5
#                 directional[i]["data"][1] = [wxy[1,1]; wxy[1,2]; -top; top]
#                 directional[i]["data"][2] = [wxy[2,1]; wxy[2,2]; -bot; dr]
#             end
#         end
#         data.data["directional"] = directional
#     end
#     return data
# end
