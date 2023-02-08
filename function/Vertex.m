%% Vertex
%  A class of functions for vertex operations.
%
%% Contribution
%  Author : Mei-Heng Yueh (yue@ntnu.edu.tw)
%  Created: 2020/02/01
% 
%  Copyright Mei-Heng Yueh
%  http://math.ntnu.edu.tw/~yueh

classdef Vertex
   methods (Static)
        function NormV = Norm(V)
            NormV = sqrt(sum(V.*V,2));
        end
        function NormV2 = Norm2(V)
            NormV2 = sum(V.*V,2);
        end
        function V = Normalize(V)
            V = bsxfun(@rdivide, V, Vertex.Norm(V));
        end
        function V = Centralize(V)
            V = bsxfun(@minus, V, mean(V, 1));
        end
        function V = Inv(V)
            V = bsxfun(@rdivide, V, Vertex.Norm2(V));
        end
        function V = Unique(V)
            V = unique(V, 'stable', 'row');
        end
        function uv = SGProj(S)
            uv = S(:,[1,2]) ./ (1-S(:,[3,3]));
        end
        function S = InvSGProj(C)
            B = 1+C(:,1).^2+C(:,2).^2;
            S = [2*C(:,1)./B, 2*C(:,2)./B, (B-2)./B];
        end
        function [InnerIdx, OuterIdx] = InnerIndex(uv, Radius)
            if nargin < 2
                Radius = 1.2;
            end
            InnerIdx = Vertex.Norm(uv) < Radius;
            if nargout == 2
                OuterIdx = ~InnerIdx;
            end
        end
   end
end