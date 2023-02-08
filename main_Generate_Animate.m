clear;clc;close all;

FileNameList = {'Horse1'; 'Horse2'; 'Horse3'; 'Horse4'; 'Horse5'; 'Horse6'; 'Horse7'; 'Horse8'};
FileNo = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1];
FileNum = length(FileNo);
M = cell(FileNum,1);
for k = 1:FileNum
    FileName = FileNameList{FileNo(k)};
    M{k} = load(fullfile('data', [FileName '.mat']));
end

F = M{1}.F;
Vno = size(M{1}.V,1);
V_Basis = zeros(Vno, 3, FileNum);

e = ones(Vno,1);
rgb = [129/255, 159/255, 247/255];
Vrgb = rgb(e,:);

% Rotation for a better visualization
theta = -pi/2;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

for k = 1:FileNum
    V_Basis(:,:,k) = M{k}.V;
    V_Basis(:,[1,3],k) = V_Basis(:,[1,3],k) * R.';
end

%%
TimeBasis = 0:FileNum-1;
dt = 5e-1;
t = 0:dt:FileNum-1;
tno = length(t);
V_Seq = zeros(Vno, 3, tno);
V_Seq(:,1,:) = spline(TimeBasis, V_Basis(:,1,:), t);
V_Seq(:,2,:) = spline(TimeBasis, V_Basis(:,2,:), t);
V_Seq(:,3,:) = spline(TimeBasis, V_Basis(:,3,:), t);

%%
camlight;
VideoName = 'Horse.avi';
VideoObject = VideoWriter(VideoName);
open(VideoObject);

PatchHandle  = patch('Faces', F, 'Vertices', V_Seq(:,:,1), 'FaceVertexCData', Vrgb, 'EdgeColor', 'none', 'FaceColor', 'interp');
ylim([-1,1]); zlim([-1,1]); xlim([-1,1]);
CurrentFrame = getframe;
writeVideo(VideoObject, CurrentFrame);

for k = 2:tno
    PatchHandle.Vertices = V_Seq(:,:,k);
    ylim([-1,1]); zlim([-1,1]); xlim([-1,1]);
    CurrentFrame = getframe;
    writeVideo(VideoObject, CurrentFrame);
end

close(VideoObject);
close all