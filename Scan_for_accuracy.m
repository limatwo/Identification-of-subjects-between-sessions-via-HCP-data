% 
clear;
% when matlab is showing an error that it cannot find the directory right
% click on the folders and select indicate files not on the path, and
% then select include folders.
% the dir command shows matlab where the data is  
%filenamesLR = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Vectorized_FC\LR/**/*.mat');
%filenamesRL = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Vectorized_FC\RL/**/*.mat');
%filenamesLR = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Data_at_1125\Vectorized_data\LR/**/*.mat');
%filenamesRL = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Data_at_1125\Vectorized_data\RL/**/*.mat');
filenamesRL = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Reduced_times_denoised_data\Data_at_417\RL/**/*.mat');
filenamesLR = dir('C:\Users\Mahmoud\OneDrive - McGill University\Bio_QLS_MCGILL\_PhD_Work\Test-Retest-reliability\Scripts\Data\Reduced_times_denoised_data\Data_at_417\LR/**/*.mat');
% mkdir ../data_at_8 LR
% 
% sourceFileName = fullfile('_Rest1_LR_ROI_data_Gordon_333_surf.mat*');
% destinationFileName = sprintf('LR');
% movefile(sourceFileName, destinationFileName);
% Also use matlab  dbstop if error for debugging
%load(filenamesLR(1).name)
%T= table (filenames);
%m=zeros(10,1);
% following script is to load data
me=0;
H=10;
Z=cell(10,1);
L=cell(10,1);
m=zeros(10);

for i = 1:length(filenamesLR)
     for j= 1:length(filenamesLR)
        holder1= load(filenamesRL(i).name);
        holder =load(filenamesLR(j).name);
        Z{i}=holder1.N;  
        L{j}=holder.N;          
     end 
     %me(i,:)=Ncorr;
end
%% To calculate correlations between LR and RL for each subject
for ii =1:10 
    for jj = 1:10
       %for tt= 1:100
       corrN=corr(Z{ii,1},L{jj,1}); % to know what is inside a loop set a breakpoint inside the loop
       %Ma_X=max(cell2mat(m));
       % b=cell2mat(corrN);
       m(ii,jj)=corrN;
       % end
         
    end
end 
% c=0;
% for k= 1:10
%     h=max(m(ii));
%     c=c+1;
%     
% end
%%
n_subj = 10;

%x = rand(n_subj,n_subj);
 
c = 0;
for i = 1:10
    tmp = m(i,:);
    [val, loc] = max(tmp);
    if loc==i
        c = c+1;
    end
end

acc = 100*c/n_subj;

% % Script for finding maximum per subject
% clear
% P=rand(10,10); 
% YY=zeros(10,1);
% 
% for kk=1:10
%     for bb=0:9
%             h=bb+1;
%             YY(bb)=(P(kk,bb)> P(kk,h));
%             if bb==9
%                 bb=1;
%             end
%            C=sum(YY);
%     end
% end
% %% 
%     A0 = 7 ;
%     B  = 1 ;
%     for k=1:10 % iteration
%       A = A0 ; % reset one variable
%       B = k*A + 2*B + k; % calculations that change another variable
%     end
% 
% % VM=load(filenamesRL(1).name)
% % L= load(filenamesLR(1).name, 'N')
% % LM= load(filenamesRL(1).name,'N')
% % C= corr(L,LM)
% % 
% % p = rand(10);
% % q = rand(10);
% % save('pqfile.mat','p','q')
% % corr(p,q);
% % corr(L,N)
% %corr('S1_LR_rest1_cor_vectorized_Mat.N','S1_RL_rest1_cor_vectorized_Mat.mat')
% % test1=corr( ,r) , load(S1_RL_rest1_cor.mat,r));
% % 
% % test1=corr( load(S1_LR_rest1_cor.mat),'r');
% % load(S1_LR_rest1_cor.mat)
% % load ( , 'data' )7\a t4load 
% % load ( files(ii).name, 'data' )
% % 
% % disp('S1_LR_rest1_cor.mat')
% %%
% % a = 10;
% % k = 0.5;
% % n = 2;
% % valueOfA = zeros(1,5);
% % for mm = 1:5
% %     a = a + (a*k) + n;
% %     valueOfA(mm) = a;
% % end