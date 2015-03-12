
EIT_fname='CR_EIT017_T_EIT_85uA_6kHz_20s_200ms_Prt1_ref1_NMBA_bw4.mat'
load(EIT_fname);
%%

 map_p=[[ 3 4 7 8 9 12 13 14 15 16 19 21 22 23 26 27 28]-1]; % Montage
   map_=[ 3 4 7 8 9 12 13 14 15 16 19 21 22 23 26 27 28]; 
    ring = [7 28	12	23	21	16	26	9	13	4	3	15	8	27	14	22];

Prt_name='Prt_1.txt';
Prt_size=15;




dZ=[];
Prt_0=[];
n=1;
for iPair=1:Prt_size
    
    for i = [1:10,12:size(EIT{iPair}.dZ,2)] 
        
        [ma,mt] = max((EIT{iPair}.dZ_avg(:,i)));
        
        if ma<200 && ...
                mt>0 && ...
                std(EIT{iPair}.dZ_avg(4500:5050,i)) <10
            
            dZ(:,n)=EIT{iPair}.dZ_avg(:,i);
           
            g=1:16;
            Prt_0(n,:) = [g(ring==EIT{iPair}.Pair(1)), g(ring==EIT{iPair}.Pair(2)), ...
                g(ring==map_(i)) ];
            n=n+1;
        end
    end
end
    
     plot(dZ);
    
 save(['for_image_' EIT_fname(1:15) '.mat'],'dZ','Prt_0','-v7.3');
    
   