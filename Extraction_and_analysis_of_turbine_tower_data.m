

%################################################################
% KROGITAL 
% Created by Likith Krishnappa
% Extraction of Turbine Tower Data, Performing Analysis and Writing it to a
% File Format that IGUS Robot can Read
% [Pefect textile]
%################################################################


clear all
clc

A=[];
B=[];
C=[];
D=[];
E=[];
F=[];
G=[];
H=[];
X=[];
Y=[];
V=[];
W=[];

  



for k= 1:1

  textFilename = ['turbine-08_2019-10-14 12_30_18+00_00_2019-10-14 15_10_43+00_00_2.csv'];    %Open text file; change the file name as needed
  fid = fopen(textFilename,'r+');

  
  A=[];                                                                     % allocate accumulator array
while ~feof(fid)
  try
    A=[A;str2num(fgetl(fid))];
    Acc_x = A(:,1);                                                         % reading accelerations
    Acc_y = A(:,2);
    Acc_z = A(:,3);
    Vel_x = A(:,4);                                                         % reading velocity                   
    Vel_y = A(:,5);
    Vel_z = A(:,6);
    Pos_x = A(:,7);                                                         % reading position
    Pos_y = A(:,8);
    Pos_z = A(:,9);
    
    
  catch
     
  end
end
end


data = 1:10000;                                                             % selecting only a part of the data 
data = data*0.033;                                                          % timestep of each data point
data = data';


%-----------------For X position--------------------------------------

Pos_x = Pos_x*1000;                                                         % converting m to mm

[f,l] = findpeaks(Pos_x,'MinPeakdistance',80);                              % find the peaks of the sine curve (upper half)
l=l*0.033;                                                                  % converting timesteps to actual time in seconds

[V,I] = findpeaks(-Pos_x,'MinPeakdistance',75);                             % find the peaks of the sine curve (lower half)
I = I*0.033;                                                                % converting timesteps to actual time in seconds

plot(data,Pos_x)                                                            % plotting the entire curve
hold on
plot(I,-V,'o')                                                              % plotting the obatined peak points
hold on
plot(l, f, 'o');                                                            % plotting the obatined peak points


f = f(1:end-1);                                                             % making sure the matrices are of same length
l = l(1:end-1);

c = [f';-V'];                                                               % combining the positive and negative peak points
c = c(:);
c = [0;c];

e = [l';I'];                                                                % combining the positive and negative time points                                                              
e = e(:);
e1 = [0;e];


%-----------------For Y position--------------------------------------


Pos_y = Pos_y*1000;                                                          % refer to the comments in the X positon section

[f1,l1] = findpeaks(Pos_y,'MinPeakDistance',70,'minpeakheight',0);
l1=l1*0.033;

[V1,I1] = findpeaks(-Pos_y,'MinPeakDistance',75,'minpeakheight',-1);
I1 = I1*0.033;

figure
plot(data,Pos_y)
hold on
plot(l1, f1, 'o');
hold on
plot(I1,-V1,'o')


V1 = V1(2:end);
I1 = I1(2:end);

c1 = [f1';-V1'];
c1 = c1(:);
c1 = [0;c1];

kl = [l1';I1'];
kl = kl(:);
kl1 = [0;kl];


%-----------------For Z position--------------------------------------

Pos_z = Pos_z*1000;                                                         % refer to the comments in the X positon section

[f2,l2] = findpeaks(Pos_z,'MinPeakDistance',70,'minpeakheight',0);
l2=l2*0.033;

[V2,I2] = findpeaks(-Pos_z,'MinPeakDistance',70,'minpeakheight',-1);
I2 = I2*0.033;

figure;
plot(data,Pos_z)
hold on
plot(I2,-V2,'o')
hold on
plot(l2, f2, 'o');


V2 = V2(2:end-1);
f2 = f2(1:end-2);

I2 = I2(2:end-1);
l2 = l2(1:end-2);

c2 = [f2';-V2'];
c2 = c2(:);
c2 = [0;c2];

gl = [l2';I2'];
gl = gl(:);
gl1 = [0;gl];
gl1 = gl1(1:end-1);


%--------------------------------------------------------------------------
%--------------------------Velocity----------------------------------------
%--------------------------------------------------------------------------


%------------------For X velocity------------------------------------------


Vel_x = Vel_x*1000;                                                          % refer to the comments in the X positon section

[p,q] = findpeaks(Vel_x,'MinPeakdistance',80,'minpeakheight',2);
q = q*0.033;

[r,s] = findpeaks(-Vel_x,'MinPeakdistance',80,'minpeakheight',-1);
s = s*0.033;

figure;
plot(data,Vel_x)
hold on
plot(q, p,'o')
hold on
plot(s, -r, 'o');

p = p(1:end-1);
q = q(1:end-1);

t = [p';r'];
t = t(:);
t = [0;t];                              % X velocity data

w = [q';s'];
w = w(:);
w = [0;w];                              % time data



%------------------For Y velocity------------------------------------------


Vel_y = Vel_y*1000;                                                          % refer to the comments in the X positon section

[p1,q1] = findpeaks(Vel_y,'MinPeakdistance',100, 'minpeakheight', 0);
q1=q1*0.033;

[r1,s1] = findpeaks(-Vel_y,'MinPeakdistance',100 ,'minpeakheight', -1);
s1 = s1*0.033;

figure
plot(data,Vel_y)
hold on
plot(q1,p1,'o')
hold on
plot(s1, -r1, 'x');


r1 = r1(2:end);
s1 = s1(2:end);

t1 = [p1';r1'];
t1 = t1(:);
t1 = [0;t1];

w1 = [q1';s1'];
w1 = w1(:);
w1 = [0;w1];                              % time data



%------------------For Z velocity------------------------------------------


Vel_z = Vel_z*1000;                                                          % refer to the comments in the X positon section

[p2,q2] = findpeaks(Vel_z,'MinPeakdistance',85,'minpeakheight',0);
q2=q2*0.033;

[r2,s2] = findpeaks(-Vel_z,'MinPeakdistance',100);
s2 = s2*0.033;

figure
plot(data,Vel_z)
hold on
plot(q2,p2,'o')
hold on
plot(s2, -r2, 'x');

p2 = p2(1:end-2);
q2 = q2(1:end-2);
r2 = r2(1:end-1);
s2 = s2(1:end-1);

t2 = [p2';r2'];
t2 = t2(:);
t2 = [0;t2];

w2 = [q2';s2'];
w2 = w2(:);
w2 = [0;w2];

resultant_vel = mean(t);                                                    % finding the mean X velocity
resultant_vel_y = mean(t1);                                                 % finding the mean Y velocity
resultant_vel_z = mean(t2);                                                 % finding the mean Z velocity

resultant_vel_2d_XY = sqrt(t.^2+t1.^2);                                     % finding the resultant XY velocity
resultant_vel_2d_XZ = sqrt(t.^2+t2.^2);                                     % finding the resultant XZ velocity
resultant_vel_2d_YZ = sqrt(t1.^2+t2.^2);                                    % finding the resultant YZ velocity


resultant_vel_3d = sqrt(t.^2+t1.^2+t2.^2);                                  % finding the resultant XYZ velocity



%--------------------------------------------------------------------------
%---------------------------acceleration----------------------------------
%--------------------------------------------------------------------------

% % Acc_x = Acc_x*1000;
% % 
% % [g,h] = findpeaks(Acc_x,'MinPeakdistance',100);
% % h=h*0.033;
% % 
% % [i,j] = findpeaks(-Acc_x,'MinPeakdistance',75);
% % j = j*0.033;
% % 
% % plot(data,Acc_x)
% % hold on
% % plot(h,g,'o')
% % hold on
% % plot(j,-i, 'o');
% % 
% % Tf = filloutliers(g,'linear');
% % [SF,TFv] = rmoutliers(i);
% % 
% % plot(h,Tf,'*')
% % hold on
% % plot(j,-SF, '*');
% % 
% % Tf= Tf(2:end);
% % Tf1 = Tf1 (2:end)
% % k = [Tf';Tf1'];
% % k = k(:);





% % Acc_y = Acc_y*1000;
% % 
% % [g1,h1] = findpeaks(Acc_y,'MinPeakdistance',80);
% % h1=h1*0.033;
% % 
% % [i1,j1] = findpeaks(-Acc_y,'MinPeakdistance',90);
% % j1 = j1*0.033;
% % 
% % plot(data,Acc_y)
% % hold on
% % plot(h1,g1,'o')
% % hold on
% % plot(j1,-i1, 'x');
% % 
% % 
% % g1= g1(1:end-2);
% % i1= i1(2:end-2);
% % k1 = [g1';i1'];
% % k1 = k1(:);
% % 
% % 
% % 
% % result_acc = sqrt(k.^2+k1.^2);


% Acc_z = Acc_z*1000;
% 
% 
% Of = isoutlier(Acc_z);
% plot(data,Acc_z)
% Of1 = filloutliers(g2,'linear');
% 
% [g2,h2] = findpeaks(Acc_z,'MinPeakdistance',80);
% h2=h2/100;
% 
% [i2,j2] = findpeaks(-Acc_z,'MinPeakdistance',80);
% j2 = j2/100
% 
% plot(data,Acc_z)
% hold on
% plot(h2,g2,'o')
% hold on
% plot(j2,-i2, 'o');
% 
% Of = filloutliers(g2,'linear');
% plot(h2,Of,'*')
% hold on
% 
% 
% g1= g1(1:end-2);
% i1= i1(2:end-2);
% k1 = [g1';i1'];
% k1 = k1(:);




%-----------------------------------------------------------------------
% -------------------  Writing Files   ------------------------------------
%---------------------------------------------------------------------------

% % % % % %----------------------------- 1D motion-----------------------------
% % % % % 
% % % % % fid=fopen('Single_speed_X.txt','w');
% % % % % fprintf(fid,'ppz z1 %f \n', [c]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % fid=fopen('Single_speed_Y.txt','w');
% % % % % fprintf(fid,'ppz z2 %f \n', [c1]');
% % % % % fclose(fid);true
% % % % % 
% % % % % fid=fopen('Single_speed_Z.txt','w');
% % % % % fprintf(fid,'ppz z3 %f \n', [c2]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % 
% % % % % fid=fopen('Variable_speed_X.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z1 %f\n\n', [t c]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % fid=fopen('Variable_speed_Y.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z2 %f \n \n', [t1 c1]');
% % % % % fclose(fid);true
% % % % % 
% % % % % fid=fopen('Variable_speed_Z.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z3 %f\n \n', [t2 c2]');
% % % % % fclose(fid);true
% % % % % 
% % % % % %------------------------------ 2D Motion---------------------------------
% % % % % 
% % % % % 
% % % % % fid=fopen('Single_speed_XY.txt','w');
% % % % % fprintf(fid,'ppz z1 %f z2 %f\n', [c c1]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % fid=fopen('Single_speed_XZ.txt','w');
% % % % % fprintf(fid,'ppz z1 %f z3 %f\n', [c c2]');
% % % % % fclose(fid);true
% % % % % 
% % % % % fid=fopen('Single_speed_YZ.txt','w');
% % % % % fprintf(fid,'ppz z2 %f z3 %f\n', [c1 c2]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % 
% % % % % fid=fopen('Variable_speed_XY.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z1 %f z2 %f\n\n', [resultant_vel_2d_XY c c1]');
% % % % % fclose(fid);true
% % % % % 
% % % % % 
% % % % % fid=fopen('Variable_speed_XZ.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z1 %f z3 %f\n \n', [resultant_vel_2d_XZ c c2]');
% % % % % fclose(fid);true
% % % % % 
% % % % % fid=fopen('Variable_speed_YZ.txt','w');
% % % % % fprintf(fid, 'ges zal kar %f \nppz z2 %f z3 %f\n \n', [resultant_vel_2d_YZ c1 c2]');
% % % % % fclose(fid);true



%--------------------------- 3D Motion------------------------------------

 c = c/2;           % halfing x values
 c1 = c1*5;         % multiplying y by 5
 c2 = c2/1.5;        

rc = diff (c);
rc1 = diff (c1);
rc2 = diff (c2);

rc = rc(1:end-1);
rc1 = rc1(1:end-1);
rc2 = rc2(1:end-1);

Nr = (3:155)';
r_resultant_vel_3d = resultant_vel_3d(2:end-1);

fid=fopen('Single_speed_XYZ.xml','w');
fprintf(fid,'<?xml version="1.0" encoding="utf-8"?> \n' )
fprintf(fid,'<Program> \n') 
fprintf(fid,'   <Header RobotName="igus Arm" RobotType="igus_5DOF_BV" GripperType="" Software="iRC V902-11-023" /> \n')
fprintf(fid,'   <Joint AbortCondition="False" Nr="1" Source="Numerical" velPercent="50" acc="40" smooth="20" a1="0" a2="0" a3="0" a4="0" a5="0" a6="0" e1="0" e2="0" e3="0" Descr="" /> \n')
fprintf(fid,'   <Relative AbortCondition="False" Nr="2" acc="40" smooth="20" MoType="cartbase" vel="100" x="0" y="0" z="100" e1="0" e2="0" e3="0" Descr="" /> \n')
fprintf(fid,'   <Relative AbortCondition="False" Nr="%d" acc="40" smooth="20" MoType="cartbase" vel="%f" x="%f" y="%f" z="%f" e1="0" e2="0" e3="0" Descr="" /> \n', [Nr r_resultant_vel_3d  rc rc1 rc2]')
fprintf(fid,'</Program>')  
fclose(fid);true






% % % % %-----------------------for noticable execution---------------------------
% % % % 
% % % % plot(c);
% % % % curve = animatedline('linewidth',2,'color','b');
% % % % time = [0:4:619]';
% % % % axis equal
% % % % set(gca, 'xlim', [-150 150],'ylim', [0 620])
% % % % for i = 1:length(c)
% % % %     addpoints(curve,c(i),time(i))
% % % %     drawnow
% % % %     pause(0.1);
% % % % end
% % % % 
% % % % fid=fopen('variable_speed.txt','w');
% % % % fprintf(fid, 'ges zal kar %f \nbes zal kar %f \nppz z1 %f z2 %f z3 %f \n \n', [resultant_vel x c c1 c2]');
% % % % fclose(fid);true
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %-----------------------------FFT of the signal--------------------------
% % % % 
% % % % NFFT=16384; %NFFT-point DFT      
% % % % X=fftshift(fft(Pos_x,NFFT)); %compute DFT using FFT      
% % % % fVals=(-NFFT/2:NFFT/2-1)/NFFT; %DFT Sample points        
% % % % plot(fVals,abs(X));      
% % % % title('Double Sided FFT - with FFTShift');       
% % % % xlabel('Normalized Frequency')       
% % % % ylabel('DFT Values');
% % % % 
% % % % fs = 1/0.033;               % number of points per second
% % % % NFFT=16384;      
% % % % X=fftshift(fft(Pos_z,NFFT));         
% % % % fVals=fs*(-NFFT/2:NFFT/2-1)/NFFT;        
% % % % plot(fVals,abs(X),'b');      
% % % % title('Double Sided FFT - with FFTShift');       
% % % % xlabel('Frequency (Hz)')         
% % % % ylabel('|DFT Values|');



% plot(Pos_x,Pos_y)
% 
% 
% write = horzcat(data,Pos_x);
% write_1 = horzcat(data,Pos_y);
% write_2 = horzcat(data,Pos_z);
% 
% write_11 = horzcat(data,Vel_x);
% write_12 = horzcat(data,Vel_y);
% write_13 = horzcat(data,Vel_z);
% 
% write_21 = horzcat(data,Acc_x);
% write_22 = horzcat(data,Acc_y);
% write_23 = horzcat(data,Acc_z);
% 
% 
% 
% filename = sprintf('pos_X_data.txt');    %to create a sequence of files
%     dlmwrite(filename,write,'delimiter', '\t')     
%     
%     
%     
% filename = sprintf('Vel_X_data.txt');    %to create a sequence of files
%     dlmwrite(filename,write_11,'delimiter', '\t')  
%     
%      xlswrite('pos_X_data.txt', write)
%      xlswrite('pos_Y_data.xlsx', write_1)
%      xlswrite('pos_Z_data.xlsx', write_2)
%      
%      myfit = fit(data,Mag_x,'sin8')
%      YHat = myfit(Mag_x)
%      
%      filename = sprintf('trial_Z_data.txt');    %to create a sequence of files
%     dlmwrite(filename,write_2,'delimiter', '\t')     




