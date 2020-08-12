%%
%filenames = dir( ' /**/*.mat')
clear
filenames = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Reduced_times_denoised_data\data_at_417/**/*.mat');
%filenames = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Data_at_1083/**/*.mat');

%load(filenames(2).name)

% H=ROI_data;
% your_mat = H(1:1125,:); 

%disp(filenames(1).name);
for i = 1:length(filenames)
    %disp(filenames(i).name);
    ROI= load(filenames(i).name);  
    Lh=(ROI.ROI_data);
    %L1= Lh;
    L1= Lh(1:417,:);
    save(filenames(i).name,'L1', '-append');
end

%% correlation calculation
for i = 1:length(filenames)
    %disp(filenames(i).name);
    G=load(filenames(i).name);
%     G=load(filenames(1).name);
%     disp(G)
    T=G.L1;
    C1=corr(T);
    save(filenames(i).name,'C1', '-append');
end
%% vectorization of the correlation matrices 
for i = 1:length(filenames)
    disp(filenames(i).name);
    Z=load(filenames(i).name);
    Z1=Z.C1;
    U=triu(Z1,1);
    J=U(:);
    a=0;
    N=setdiff(J,a,'stable');
    save(filenames(i).name,'N', '-append') ;
end
%% PLOT 
clear
% denoised 
x=[0.5;0.8;1;2;3;5;7;9;11;13;13.5;14.4];
y=[10;50;60;70;70;80;70;80;70;80;80;80];
plot(x,y)
xlabel('Time in minutes') 
ylabel('Accuracy in percentage') 
title('with denoised data')
clear
z=[0.5;1;2;3;5;7;9;11;13;13.5;14.4];
t=[20;30;30;40;50;60;60;60;70;70;70];
plot(z,t)
xlabel('Time in minutes') 
ylabel('Accuracy in percentage') 
title('with non processed data')
%%% Load data from two different files
% path_root = pwd; 
% 
% for s = 1 : 10
%       path_subj = [path_root,'\Data_at_1125',num2str(s),'_LR_ROI_data_Gordon_333_surf_Mat'];
%      %100206_Rest1_LR_ROI_data_Gordon_333_surf
%      load(path_subj);
% end 
% 
% %% Correlation Calculation

for i = 1:length(filenames)
    N1=load(filenames(i).name);
    N2=N1.N;
end







