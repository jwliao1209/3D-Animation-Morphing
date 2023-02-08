function Vnew = MeshEditing(F, V, landmark, V_landmark)
Vno  = size(V,1);
L    = Tri.Laplacian(F,V);
Lno  = length(landmark);
C    = sparse((1:Lno).', landmark, 1, Lno, Vno);
A    = [L; C];
rhs  = [L*V; V_landmark];
Vnew = A\rhs;