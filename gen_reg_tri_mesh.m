function [tri, nodnum, xcoord, ycoord, e2d, n2d] =  gen_reg_tri_mesh(x, y)

%%%    Given array in x and y, generates triangular mesh 
%%%    with x,y coordinates, number of points, and distances
%%%    for each triangle edge. 

lon = x;
lat = y;

ddx = x(2) - x(1);%x(2) - x(1);
ddy = y(2) - y(1);%y(2) - y(1);

%%%    Borrowed from Sergey Donilov, AWI May 2025
            
% Channel: regular mesh, equilateral triangles
     
  %ddy=ddx*sqrt(3)/2;

  %lon=6400000*pi*(0:ddx:10)/180;
  %lat=6400000*pi*(0:ddy:10)/180;
  
  % % distorted mesh
  % if distorted==1
  %     lon = zeros(1);
  %     k = 0;
  %     dx=0;
  %     while k < 10
  %         dx=ddx + 0.01*sin(10*pi*6400000*pi*k/180);
  %         k=k+dx;
  %         lon=[lon; [k]];
  %     end
  %     lon = (6400000*pi*lon/180)';
  % end
  
  
  
  
  
  nx=length(lon);
  ny=length(lat);

  nodnum=1:nx*ny;
  nodnum=reshape(nodnum,[ny, nx]);

  xcoord=zeros([ny, nx]);
  xnum=zeros([ny, nx]);

  ycoord=xcoord;
  ynum=zeros([ny, nx]);

  for n=1:nx
      ycoord(:,n)=lat';
      ynum(:,n)=(1:ny)';
  end
   
  for n=1:ny
      xcoord(n,:)=lon;
      xnum(n,:)=(1:nx);
  end

  for n=2:2:ny
    xcoord(n,:)=xcoord(n,:)+0.5*ddx;
  end

  xcoord=reshape(xcoord,[1,nx*ny]);
  ycoord=reshape(ycoord,[1,nx*ny]); 
  
  tri=[];
  for n=1:nx-1,
      for nn=1:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
      tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
      end;
      for nn=2:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn+1,n+1)]];
      tri=[tri; [nodnum(nn,n),nodnum(nn,n+1),nodnum(nn+1,n+1)]];
      end
  end


v_grid = NaN(nx,ny);

counter = 0;
for j = 1:ny
    for i = 1:nx
        counter = counter + 1;
        v_grid(i,j) = counter;
    end
end

  tri=tri';
  e2d=length(tri(1,:));
  n2d=nx*ny;
  
  
  % derivatives:
  dx=zeros([3, e2d]);
  dy=zeros([3, e2d]);
  
  for n=1:e2d
      x2=xcoord(tri(2,n))-xcoord(tri(1,n));
      x3=xcoord(tri(3,n))-xcoord(tri(1,n));
      y2=ycoord(tri(2,n))-ycoord(tri(1,n));
      y3=ycoord(tri(3,n))-ycoord(tri(1,n));
      d=x2*y3-y2*x3;
      dx(1,n)=(-y3+y2)/d;
      dx(2,n)=y3/d;
      dx(3,n)=-y2/d;
      
      dy(1,n)=-(-x3+x2)/d;
      dy(2,n)=-x3/d;
      dy(3,n)=x2/d;
  end;     


  dx;
      
