function [edges,edge_num,elem_edges,edge_tri,ed2d,ed2d_in,edge_dxdy, edge_cross_dxdy] = get_fv_arrays(tri, xcoord, ycoord, n_num,ng_num, nghbr_n2d, n2d,nod_in_e2d, e2d)

% Auxiliary arrays for finite-volume advection methods

% We have:  n_num
%           nod_in_e2d
%           ng_num
%           nghbr_n2d
%  These are arrays of neighboring elements (triangles) and vertices
%  (nodes) 

% Find edges and triangles containing them.
%

edge_num=0;

  for n=1:n2d
      % form edges with n by cycling over neighbor
      % nodes (if edges are not accounted yet). 
      % New edges are added only if neighbor>n  
	  
   for q=1:ng_num(n)
       node=nghbr_n2d(q,n);
	  if node>n, 
	  edge_num=edge_num+1;   % new edge (n,node)
          end;	  
   end;
   end;


   edges=zeros([2,edge_num]);
   edge_tri=zeros([2,edge_num]);
   edge_num=0;
   % Repeat and fill:
   for n=1:n2d
     for q=1:ng_num(n)
       node=nghbr_n2d(q,n);
	  if node>n, 
	  edge_num=edge_num+1;   % new edge (n,node)
	  edges(1,edge_num)=n;
	  edges(2,edge_num)=node;
	  % find triangles containing n and node
	  flag=0;
	  for k=1:n_num(n),
	  elem=nod_in_e2d(k,n);
	  elnodes=tri(:,elem);
	    for q1=1:3,
		   if elnodes(q1)==node,
		   flag=flag+1;
		   edge_tri(flag,edge_num)=elem;
		   end
	    end
	  end
	  % flag=1 if second triangle is absent (boundary edge)
	  if (flag>2)|(flag==0), disp('flag'); end;
	  
	  end
      end
   end
   

   % Sort edges so that first we have those with two
   % element-neighbors.
   ai=find(edge_tri(2,:)>0);
   ed2d_in=length(ai);
   aux1=edge_tri(:,ai);
   aux2=edges(:,ai);
   ai=find(edge_tri(2,:)==0);
   aux3=edge_tri(:,ai);
   aux4=edges(:,ai);
   ed2d=ed2d_in+length(ai);
   edge_tri=[aux1,aux3];
   edges=[aux2,aux4];
   clear aux1 aux2 aux3 aux4

   % Edges from ed2d_in+1 to ed2d lie on the horizontal boundary
   % The rest (1:ed2d_in) are internal edges


   % the list of elements on both sides of edge e, edge_tri(:,e)
   % should be ordered so that the first is a triangle which is to the left of 
   % the edge (vector pointing from the first to the second node of the edge. 
   % If the edge is on the boundary, there is only the first triangle.  

   for n=1:ed2d
   ed=edges(:,n);
   xc(1)=sum(xcoord(tri(:,edge_tri(1,n))))/3.0-xcoord(ed(1));
   xc(2)=sum(ycoord(tri(:,edge_tri(1,n))))/3.0-ycoord(ed(1));
   xe(1)=xcoord(ed(2))-xcoord(ed(1));
   xe(2)=ycoord(ed(2))-ycoord(ed(1));
   
     if xc(1)*xe(2)-xc(2)*xe(1)>0,
	 % Vector drawn to the center of the first triangle is to the right
	 % of the edge vector. Triangles have to be exchanged:
	 elem=edge_tri(1,n);
	 elem1=edge_tri(2,n);
	 if elem1>0,
	   edge_tri(1,n)=elem1;
	   edge_tri(2,n)=elem;
	 else
	   elem=edges(2,n);         % change the order of nodes
	   edges(2,n)=edges(1,n);
	   edges(1,n)=elem;
   	 end
     end
   end

% We need an array inverse to edge_tri listing edges
% of a given triangle 
elem_edges=zeros([3,e2d]);
aux=zeros([1,e2d]);

for n=1: ed2d,
   for k=1:2,
      q=edge_tri(k,n);   % triangle number
	  if (q>0),
	  aux(q)=aux(q)+1;
	  elem_edges(aux(q),q)=n;
	  end
   end
end
 clear aux
% The edges in this list should be ordered so that they
% are listed in the same rotation sense as nodes.
for elem=1:e2d,
   elnodes=tri(:,elem);
   eledges=elem_edges(:,elem);
   for q=1:3
      for k=1:3
         if((edges(1,eledges(k))~=elnodes(q)) & ...
            (edges(2,eledges(k))~=elnodes(q))), 
           elem_edges(q,elem)=eledges(k);
	   
         end 
      end
   end
end


% The edge and elem lists agree in the sense that edge1 does not
% contain node 1 and so on

% ============================================================
% We have got:  edges(2,ed2d)
%                edge_tri(2,ed2d)
%		elem_edges(3,e2d)
%		ed2d_in, ed2d 
% ============================================================

 edge_dxdy=zeros([2,ed2d]);
 edge_cross_dxdy=zeros([4,ed2d]);
  
 % ===========
 % Distances along the edge
 % ===========
 for n=1:ed2d,
 ed=edges(:,n);
 a(1)=xcoord(ed(2))-xcoord(ed(1));
 a(2)=ycoord(ed(2))-ycoord(ed(1));
 edge_dxdy(:,n)=a;
 end;

 % ===========
 % Cross-distances for the edge
 % ===========
 for n=1: ed2d,
 ed=edges(:,n);
 el=edge_tri(:,n);
 b(1)=sum(xcoord(tri(:,el(1))))/3.0-sum(xcoord(ed))/2.0;
 b(2)=sum(ycoord(tri(:,el(1))))/3.0-sum(ycoord(ed))/2.0;
 edge_cross_dxdy(1:2,n)=b;
 if(el(2)>0),
 b(1)=sum(xcoord(tri(:,el(2))))/3.0-sum(xcoord(ed))/2.0;
 b(2)=sum(ycoord(tri(:,el(2))))/3.0-sum(ycoord(ed))/2.0;
 edge_cross_dxdy(3:4,n)=b;
 else
 edge_cross_dxdy(3:4,n)=0.0;
 end
 end


% % up/dn edge triangles (for muscl-type advection scheme)
% 
% edge_up_dn_tri=zeros([2,ed2d]);
% x=zeros([1,2]);
% b=x;
% x=x;
% for n=1: ed2d,
%    ednodes=edges(:,n); 
%    x(1)=xcoord(ednodes(2))-xcoord(ednodes(1));
%    x(2)=ycoord(ednodes(2))-ycoord(ednodes(1));
%    % Find upwind (in the sense of x) triangle, i. e. 
%    % find which triangle contains -x:
%    x=-x;
%    for k=1:n_num(ednodes(1)),
%       elem=nod_in_e2d(k,ednodes(1));
%       elnodes=tri(:,elem);
%       if (elnodes(1)==ednodes(1)), 
%          b(1)=xcoord(elnodes(2))-xcoord(elnodes(1));
% 	 b(2)=ycoord(elnodes(2))-ycoord(elnodes(1));
% 
% 	 c(1)=xcoord(elnodes(3))-xcoord(elnodes(1));
%  	 c(2)=ycoord(elnodes(3))-ycoord(elnodes(1));
%       elseif (elnodes(2)==ednodes(1)),
%          b(1)=xcoord(elnodes(1))-xcoord(elnodes(2));
% 	 b(2)=ycoord(elnodes(1))-ycoord(elnodes(2));
% 
% 	 c(1)=xcoord(elnodes(3))-xcoord(elnodes(2));
%  	 c(2)=ycoord(elnodes(3))-ycoord(elnodes(2));
%       else	 
%          b(1)=xcoord(elnodes(1))-xcoord(elnodes(3));
% 	 b(2)=ycoord(elnodes(1))-ycoord(elnodes(3));
% 
% 	 c(1)=xcoord(elnodes(2))-xcoord(elnodes(3));
%  	 c(2)=ycoord(elnodes(2))-ycoord(elnodes(3));
%       end;
%       % the vector x has to be between b and c
%       % Decompose b and x into parts along c and along (-cy,cx), i.e.
%       % 90 degree counterclockwise
%       cr=sum(c.*c);
%       bx=sum(b.*c)/cr;
%       by=(-b(1)*c(2)+b(2)*c(1))/cr;
%       xx=sum(x.*c)/cr;
%       xy=(-x(1)*c(2)+x(2)*c(1))/cr;
%       ab=atan2(by,bx);
%       ax=atan2(xy,xx);
%       % Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
%       % be what is needed
%       if (ab>0) & (ax>0) & (ax<ab),
%       edge_up_dn_tri(1,n)=elem;
%       end; 
%       if (ab<0) & (ax<0) & (ax>ab),
%       edge_up_dn_tri(1,n)=elem;
%       end;
%       if (ab==ax) | (ax==0.0),
%       edge_up_dn_tri(1,n)=elem;
%       end;
%    end;
%    % Find downwind element
%    x=-x;
%    for k=1:n_num(ednodes(2)),
%       elem=nod_in_e2d(k,ednodes(2));
%       elnodes=tri(:,elem);
%       if elnodes(1)==ednodes(2),
%          b(1)=xcoord(elnodes(2))-xcoord(elnodes(1));
% 	 b(2)=ycoord(elnodes(2))-ycoord(elnodes(1));
% 	 c(1)=xcoord(elnodes(3))-xcoord(elnodes(1));
% 	 c(2)=ycoord(elnodes(3))-ycoord(elnodes(1));
%       elseif elnodes(2)==ednodes(2)
%          b(1)=xcoord(elnodes(1))-xcoord(elnodes(2));
% 	 b(2)=ycoord(elnodes(1))-ycoord(elnodes(2));
% 	 c(1)=xcoord(elnodes(3))-xcoord(elnodes(2));
% 	 c(2)=ycoord(elnodes(3))-ycoord(elnodes(2));
%       else	 
%          b(1)=xcoord(elnodes(1))-xcoord(elnodes(3));
% 	 b(2)=ycoord(elnodes(1))-ycoord(elnodes(3));
% 	 c(1)=xcoord(elnodes(2))-xcoord(elnodes(3));
% 	 c(2)=ycoord(elnodes(2))-ycoord(elnodes(3));
%       end;
%       % the vector x has to be between b and c
%       % Decompose b and x into parts along c and along (-cy,cx), i.e.
%       % 90 degree counterclockwise
%       cr=sum(c.*c);
%       bx=sum(b.*c)/cr;
%       by=(-b(1)*c(2)+b(2)*c(1))/cr;
%       xx=sum(x.*c)/cr;
%       xy=(-x(1)*c(2)+x(2)*c(1))/cr;
%       ab=atan2(by,bx);
%       ax=atan2(xy,xx);
%       % Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
%       % be what is needed
%       if (ab>0) & (ax>0) & (ax<ab),
%       edge_up_dn_tri(2,n)=elem;
%       end;
%       if (ab<0) & (ax<0) & (ax>ab),
%       edge_up_dn_tri(2,n)=elem;
%       end;
%       if (ab==ax) | (ax==0.0),
%       edge_up_dn_tri(2,n)=elem;
%       end;
%    end;
% end;
% 
% 
% % There is problem with edges close to boundary --- there can be no  
% % up or downwind elements. We will return to the standard linear
% % upwind at nodes that belong to such edges. 
% 
% % Count the number of 'good' edges:
% k=0; 
% for n=1:ed2d
%    if (edge_up_dn_tri(1,n)~=0) & (edge_up_dn_tri(2,n)~=0), k=k+1;
% end
% end;
% disp('UP-DN ')
%  k
%  ed2d   
% 
% % Higher-order MUSCL: on up/dn elements find nodes of
% % opposing edge and coefficients of reconstruction:
% edge_up_dn_nodes=zeros([4,ed2d]);
% edge_up_dn_coef=zeros([2,ed2d]);
% for n=1:ed2d
%     ednodes=edges(:,n);
%     for j=1:2,
%     eup=edge_up_dn_tri(j,n);
%     if eup~=0, 
%        elnodes=tri(:,eup);
%        ai=find(elnodes~=ednodes(j));
%        edge_up_dn_nodes(1+(j-1)*2:2+(j-1)*2,n)=elnodes(ai);
%        % find the intersection between the edge n and the edge containing 
%        % the nodes (it should be one of these nodes on uniform meshes)
%         x(1)=xcoord(ednodes(2))-xcoord(ednodes(1));
%         x(2)=ycoord(ednodes(2))-ycoord(ednodes(1)); 
%        % x is the along the edge 
%         b(1)=xcoord(elnodes(ai(2)))-xcoord(elnodes(ai(1)));
%         b(2)=ycoord(elnodes(ai(2)))-ycoord(elnodes(ai(1)));
%        % b is the vector to be intersected 
%        % find their intersection and interpolation coefficients
%        %%% A new implementation 
%         r1=xcoord(ednodes(1))-xcoord(elnodes(ai(1)));
%         r2=ycoord(ednodes(1))-ycoord(elnodes(ai(1)));
%         r1=r1*x(2)-r2*x(1);        
%         r2=b(1)*x(2)-b(2)*x(1);
%         r1=r1/r2;
%         edge_up_dn_coef(j,n)=max(0.0,min(1.0,1.0-r1));
%        %%%%% END of new
%        if 0>1, % old piece
%         r1=ycoord(ednodes(1))*x(1)-xcoord(ednodes(1))*x(2);
%         r2=ycoord(elnodes(ai(1)))*b(1)-xcoord(elnodes(ai(1)))*b(2);
%         dd=x(1)*b(2)-x(2)*b(1);
%         x_int=(r1*b(1)-r2*x(1))/dd;
%         y_int=(r1*b(2)-r2*x(2))/dd;
%         if abs(b(1))> 0.0001 
%         c_int=abs((xcoord(elnodes(ai(2)))-x_int)/b(1));
%         edge_up_dn_coef(j,n)=c_int;
%         else
%         c_int=abs((ycoord(elnodes(ai(2)))-y_int)/b(2));
%         edge_up_dn_coef(j,n)=c_int;
%         end;
%        end   % end old piece
%     end;
%     end;
% end;
% % cint is the coefficient for the contribution from the first node
% % edge_up_dn_nodes(1...) for up-edge and edge_up_dn_nodes(3...) for dn 
% % max(max(edge_up_dn_coef))
% % min(min(edge_up_dn_coef))
% 
% % ==========================
% % Derivatives of nc basis functions:
% % ==========================
% bafuNCx=zeros([3,e2d]);
% bafuNCy=zeros([3,e2d]); 
% for elem=1:e2d 
%    ed=elem_edges(:,elem);
%    elnodes=tri(:,elem);
%    for k=1:3
%       a=[xcoord(edges(1,ed(k))), ycoord(edges(1,ed(k)))];
%       b=[xcoord(edges(2,ed(k))), ycoord(edges(2,ed(k)))];
%       b=b-a;
%             % The derivative for triangle elem and ed(k) is
%       % -2/h in the direction of inner normal drawn from ed(k) 
%       n=[-b(2), b(1)]; % inner normal for the triangle on the left:
%       if edge_tri(2,ed(k))==elem, n=-n; end; 
%       bafuNCx(k,elem)=-n(1)/elem_area(elem);
%       bafuNCy(k,elem)=-n(2)/elem_area(elem);
%    end
% end
% 
%  % derivatives of P1 basis functions:
%   dx=zeros([3, e2d]);
%   dy=zeros([3, e2d]);
% 
%   for n=1:e2d,
%       x2=xcoord(tri(2,n))-xcoord(tri(1,n));
%       x3=xcoord(tri(3,n))-xcoord(tri(1,n));
%       y2=ycoord(tri(2,n))-ycoord(tri(1,n));
%       y3=ycoord(tri(3,n))-ycoord(tri(1,n));
%       d=x2*y3-y2*x3;
%       dx(1,n)=(-y3+y2)/d;
%       dx(2,n)=y3/d;
%       dx(3,n)=-y2/d;
% 
%       dy(1,n)=-(-x3+x2)/d;
%       dy(2,n)=-x3/d;
%       dy(3,n)=x2/d;
%   end;     
%   % =====================
%   % cell neighbors and cell gradients 
%   % =====================
%   % =============
%   % cell neighbors == the cells that share edges
%   % =============
%    ee_pos=zeros([3,e2d]);
%    for elem=1:e2d,
%    eledges=elem_edges(:,elem);
%    for j=1:3,
%    elem1=edge_tri(1,eledges(j));
%    if elem1==elem, 
%        elem1=edge_tri(2,eledges(j)); 
%        if elem1<1, elem1=elem; end
%    end
%    ee_pos(j,elem)=elem1;
%    end; end
%   % =============== 
%   % cell gradients
%   % ===============
%    rc=[0,0]; xn=[0,0,0]; yn=xn;
%    g_vec=zeros([6,e2d]);
%    for elem=1:e2d,
%      rc(1)=sum(xcoord(tri(:,elem)))/3;
%      rc(2)=sum(ycoord(tri(:,elem)))/3;
%      for j=1:3,
%         elemn=ee_pos(j,elem);
%         if elemn>0, 
%            xn(j)=sum(xcoord(tri(:,elemn)))/3-rc(1);
%            yn(j)=sum(ycoord(tri(:,elemn)))/3-rc(2);
%         end
%      end
%      cxx=sum(xn.^2);
%      cxy=sum(xn.*yn);
%      cyy=sum(yn.^2);
%      dfactor=cxy*cxy-cxx*cyy;
%      if abs(dfactor)>0.001*(cxx+cyy),
%      g_vec(1:3,elem)=(cxy*yn-cyy*xn)/dfactor;
%      g_vec(4:6,elem)=(cxy*xn-cxx*yn)/dfactor;
%      end
%    end
