function problem2
clear all; close all;
nc = 45; %number of fixed points on the circle
a1 = 10;
tol = 10^-5;
phi = (1:nc)'*(2*pi/nc);
cfix = 1.5+[cos(phi),sin(phi)];
fd = @(p) drectangle(p,0,3,0,3);
[pts,tri] = distmesh2d(fd,@huniform,0.075,[0,0;3,3],[0,0;0,3;3,0;3,3;cfix]);

D0 = find(pts(:,1)<tol);
D1 = find(pts(:,1)>3-tol);
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,union(D0,D1));
A = sparse(Npts,Npts);
b = sparse(Npts,1);

for j =1:Ntri
    if in(pts(tri(j,1)),pts(tri(j,2)),pts(tri(j,3)))       
        A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:))+ a1*stima3(pts(tri(j,:),:));
    else        
        A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:))+ stima3(pts(tri(j,:),:));
    end
    b(tri(j,:))=0;
end
u=sparse(Npts,1);
u(D1)=1;
b=b-A*u;
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

abs_current_centers = zeros(Ntri,1);
for j = 1:Ntri
    vertices = pts(tri(j,:),:);
    H = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
    nablau = [u(tri(j,1)),u(tri(j,2)),u(tri(j,3))]*H;
    if in(pts(tri(j,1)),pts(tri(j,2)),pts(tri(j,3)))
        abs_current_centers(j)=(a1+0.5)*norm(nablau);
    else        
        abs_current_centers(j)=norm(nablau);
    end
end

abs_current_verts = zeros(Npts,1);
count_tri = zeros(Npts,1);
for j = 1:Ntri
abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) + abs_current_centers(j);
count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
end
abs_current_verts = abs_current_verts./count_tri;

figure;
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
hold on
axis ij
view(2)

figure;
trisurf(tri,pts(:,1),pts(:,2),full(abs_current_verts)','facecolor','interp')
hold on
axis ij
view(2)

function M = stima3(vertices)
    d = size(vertices,2);
    G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
    M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
    if in(vertices(1,:),vertices(2,:),vertices(3,:))        
        M = a1*M;
    end
end

function flag = in(a,b,c)    
    %if one vertex is in the circle then the triangle should be in
    %the circle.
    if (norm(a-1.5)<1-tol|norm(b-1.5)<1-tol)|norm(c-1.5)<1-tol
        flag = true;
    else
        flag = false;
    end
end
end