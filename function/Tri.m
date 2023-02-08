%% Tri
%  A class of functions for triangular mesh operations.
%
%% Contribution
%  Author : Mei-Heng Yueh (yue@ntnu.edu.tw)
%  Created: 2020/02/01
% 
%  Copyright Mei-Heng Yueh
%  http://math.ntnu.edu.tw/~yueh

classdef Tri
  methods (Static)
    function [E12, E23, E31] = HalfEdge(F, V)
      E12 = V(F(:,2),:) - V(F(:,1),:);
      E23 = V(F(:,3),:) - V(F(:,2),:);
      E31 = V(F(:,1),:) - V(F(:,3),:);
    end
      
    function NF = Normal(F, V)
      [Vno, Dim] = size(V);
      if Dim == 2
        V = [V, zeros(Vno,1)];
      end
      E12 = V(F(:,2),:) - V(F(:,1),:);
      E13 = V(F(:,3),:) - V(F(:,1),:);
      NF = cross(E12, E13);
      NF = Vertex.Normalize(NF);
    end
    
    function A = Angle(V, F)
      [Vno, Dim] = size(V);
      if Dim == 2
        V = [V, zeros(Vno,1)];
      end
      E1 = V(F(:,2),:)-V(F(:,3),:);
      E2 = V(F(:,3),:)-V(F(:,1),:);
      E3 = V(F(:,1),:)-V(F(:,2),:);
      E1 = Vertex.Norm(E1);
      E2 = Vertex.Norm(E2);
      E3 = Vertex.Norm(E3);
      Fno = size(F,1);
      A = zeros(Fno,3);
      A(:,1) = acos( ( E2.^2 + E3.^2 - E1.^2 ) ./ ( 2.*E2.*E3 ) );
      A(:,2) = acos( ( E1.^2 + E3.^2 - E2.^2 ) ./ ( 2.*E1.*E3 ) );
      A(:,3) = acos( ( E1.^2 + E2.^2 - E3.^2 ) ./ ( 2.*E1.*E2 ) );
    end
    
    function AD = AngleDiff(F, V, U)
      AngleV = Tri.Angle(V, F);
      AngleU = Tri.Angle(U, F);
      AngleDiff = abs(AngleV - AngleU);
      AD = rad2deg(AngleDiff);
    end
    
    function A = Area(F, V)
      if size(V,2) == 2
        V = [V, 0*V(:,1)];
      end
      V12 = V(F(:,2),:) - V(F(:,1),:);
      V13 = V(F(:,3),:) - V(F(:,1),:);
      Z = cross(V12, V13);
      A = 0.5*Vertex.Norm(Z);
    end
    
    function V = AreaNormalize(F, V)
      V = Vertex.Centralize(V);
      A = Tri.Area(F, V);
      V = V ./ sqrt(sum(A));
    end
    
    function [VB, VI] = Boundary(F)
      Vno = max(max(F));
      if nargin == 1
        V = zeros(Vno,2);
      end
      M = triangulation(F, V);
      VB = freeBoundary(M);
      VB = VB(:,1);
      if nargout > 1
        VI = setdiff((1:Vno).', VB);
      end
    end
    
    function P = Plot(F, V, U)
      NF = Tri.Normal(F,V);
      Vno = size(V,1);
      e = ones(Vno,1);
      rgb = [129/255, 159/255, 247/255];
      Vrgb = rgb(e,:);
      if exist('U', 'var')
        P = patch('Faces', F, 'Vertices', U, 'FaceVertexCData', Vrgb, 'EdgeColor','w','FaceColor','interp', 'EdgeAlpha', 0.5, 'EdgeLighting', 'flat');
      else
        P = patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Vrgb, 'EdgeColor','none','FaceColor','interp', 'EdgeAlpha', 0.5, 'EdgeLighting', 'flat');
        P.FaceLighting = 'phong';
      end
      P.FaceNormals = -NF;
      camlight('headlight');
      light('Position', [1,1,1]);
      light('Position', -[1,1,1]);
      set(gcf, 'color', [0 0 0]);
      axis equal off
    end
    
    function Surf(F, V, Vrgb)
      if ~exist('Vrgb', 'var')
        Vno = size(V,1);
        e = ones(Vno,1);
        rgb = [129/255, 159/255, 247/255];
        Vrgb = rgb(e,:);
        camlight;
      end
      patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Vrgb, 'EdgeColor', 'none', 'FaceColor', 'interp');
      axis equal off
    end
    
    function Q = Quality(F, V)
      [E12, E23, E31] = Tri.HalfEdge(F, V);
      LE12 = Vertex.Norm(E12);
      LE23 = Vertex.Norm(E23);
      LE31 = Vertex.Norm(E31);
      E = [LE12, LE23, LE31];
      Q = Vertex.Norm(bsxfun(@minus, E, mean(E, 2)));
    end
    
    function [Fid, Vid] = Center(F, V)
      L = Tri.Laplacian(F, V);
      Vno = size(V,1);
      I = speye(Vno);
      dt = 1e4;
      A = I + dt*L;
      VB = Tri.Boundary(F);
      f = zeros(Vno,1);
      f(VB) = 1;
      f = A\f;
      [~, Vid] = min(f);
      Q = Tri.Quality(F, V);
      Fid = find(sum(ismember(F, Vid),2)>0);
      [~, Qid] = min(Q(Fid));
      Fid = Fid(Qid);
    end
    
    function [L, K] = Laplacian(F, V, Sigma)
      if ~exist('Sigma', 'var')
          Sigma = 1;
      end
      Vno = size(V,1);
      if size(V,2)==2
          V = [V, zeros(Vno,1)];
      end
      [E12, E23, E31] = Tri.HalfEdge(F, V);
      Fno = size(F,1);
      W = zeros(Fno,3);
      W(:,1) = -0.5*sum(E31.*E23, 2) ./ sqrt( sum(cross(E31,E23).^2, 2) ) ./ Sigma;
      W(:,2) = -0.5*sum(E12.*E31, 2) ./ sqrt( sum(cross(E12,E31).^2, 2) ) ./ Sigma;
      W(:,3) = -0.5*sum(E23.*E12, 2) ./ sqrt( sum(cross(E23,E12).^2, 2) ) ./ Sigma;
      K = sparse(F, F(:,[2, 3, 1]), W, Vno, Vno);
      K = K + K.';
      L = diag( sum(K, 2) ) - K;
    end
    
    function E = Energy(F, V, L)
      if ~exist('L', 'var')
        L = Tri.Laplacian(F, V);
      end
      Area = Tri.Area(F, V);
      H = L*V;
      E = 0.5*sum(sum(V.*H));
      E = E - sum(Area);
    end
  end
end