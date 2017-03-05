clear all
close all
clc
SNRdB=1:1:12;                       %SNR in dB
SNR=10.^(SNRdB/10);
ol=[2 1 3 1];
y1=de2bi(ol);                      %decimal to binary conversion
y1=fliplr(y1);
[r1 c1]=size(y1);
y2=reshape(y1',1,r1*c1);

[r2 c2]=size(y2);
l=mod(c2,3);
l1=3-l;                            
y3=[y2 zeros(1,l1)];
l2=length(y3);
y4=zeros(length(y3)/3,3);          % Grouping
[r c]=size(y4);
l3=l2/3;
k=1;
for i=1:r
    for j=1:c
    y4(i,j)=y3(k);
    k=k+1;
    end
end



for i=1:length(y3)
    if y4(i)==0
        y4(i)=-1;
    else 
        y4(i)=1;
    end
end
y6=zeros(r,16);                          %making_signal_periodic
 [k x]=size(y6);
%0 for i=1:k
%   
%       for j=1:16
%           y6(i,j)=y5(1,i);
%   
%       end
% end
y6=reshape(y6',k*x,1);
y6=y6';
 a1 = [-1 1 1 1 1 -1 1 -1 1 1 -1 1 1 1 -1 1 1 1 1 -1 1 1 1 1 1 1];
% pulseshaper=rcosfir(0.5, 1, 16, 1);
% [w c]=grpdelay(pulseshaper);
% w=mean(w);
% for i=1:length(a)
%     
%     if (a(i)==0);
%         a(i)=1;
%     end
%     
% end
a6=1;
for z=1:length(a1)                       %upsampling
  
      for j=1:16
          if(j==1)
          upsampled(a6)=a1(z);
          else
          upsampled(a6)=0;
          
          end
         a6=a6+1;
          
      end
end
[x2 x3]=size(upsampled);
upsampled=reshape(upsampled',1,x2*x3);
upsampling_size=16*length(a1);
 upsampling_size1=upsampling_size+16;
% a6=1;a8=1;
%      for a4=1:upsampling_size1
%     
%         if (a6<17)
%         
%             upsampled_bits1(a6)=0;
%         
%         else
%         
%             upsampled_bits1(a6)=upsampled(a8);
%             a8=a8+1;
%         end
%         a6=a6+1;
%     
%         
%      end

%      a6=1;a8=1;
%      for a4=1:upsampling_size1
%     
%         if (a6>upsampling_size)
%         
%             upsampled_bits2(a6)=0;
%         
%         else
%         
%             upsampled_bits2(a6)=upsampled(a8);
%             a8=a8+1;
%         end
%         a6=a6+1;
%     
%         
%      end
%      
% overlapped=upsampled_bits2+upsampled_bits1;
%convolution of pulse shaper and periodic signal
% [s1 z1]=size(upsampled);
% y101=reshape(upsampled',1,s1*z1);
% y101=y101*4*pi*(.25);
% y101=[y101 0];



% x2=zeros(1,79);
% for i=1:79
%     x2(i+1)=i;
%     if(i>16)
%         x2(i+1)=16;
%     end
% end


% x1=length(x);
% k=2;
% k3=2;
% finalphase(1) = 0;

k2=1;
k3=2;
finalphase(1)=0;
for i1=1:upsampling_size              %Phase Generator
  
      
  
     
      finalphase(k3)= upsampled(i1) + finalphase(k3-1);
      k3=k3+1;
%         acc = acc + val;
end
for i5=1:upsampling_size
    
    finalphase1(1,i5)=finalphase(1,i5+1);
end

 yt11=(finalphase1*2*pi*(.5))/2;             %2 *pi * h where h is modulation index with a value 1/2

for i1=1:length(yt11)
    yt11(i1)=(yt11(i1)*57.2958);                          
    if yt11(i1)>=360
        yt11(i1)=(yt11(i1)-360);
    end
    
end

clear i
for k=1:length(yt11)                        %imaginary and real signal to be transmitted
I(k)=cosd(yt11(k));
Q(k)=sind(yt11(k));
transmitted1(k)=I(k)+i*Q(k);
end



phase=[0,90,180,270];

 
 
for i=1:length(phase)                               %MLSE reciever  
  xx(1,i)=cosd(phase(1,i));
  yy(1,i)=sind(phase(1,i));
     
    
end
clear i;
phase1=xx+i*yy;
phase2=conj(phase1);                           % conjugate of phase

b=0;
b1=0;
kk=1;
for i=1:length(phase2)
    for j=1:length(transmitted1)
        a=phase2(i)*transmitted1(j);
        
        b=b+a;
        
        if(mod(j,16)==0)
            arr(kk)=b;
            
            kk=kk+1;
            b=0;
            
        end
    end
    cell{i}=arr;
    
   kk=1;
   ui=(cell2mat(cell))/16;
      
end
length_1=length(ui)/4;

% for i=1:4
%     
%         cell{i}=cell{i}/16;
%         
%    
% end

kkk=1;
km=1;
d=1;
p=0;
for i9=1:length_1
   
  
  for j2=1:4
     bbn(1,j2)=ui(1,d+p);
       
  arr12(i9,j2)=bbn(1,j2);
  vv(1,j2)=arr12(i9,j2);
  d=d+length_1;
  if (j2==4)
%    ard(i9)=find(max(real(vv)));
   [actt(i9),actt1(i9)]=max(real(vv));                    
   
  end
  
  
  end
   
   
  d=1;
  p=p+1;
end




    

%     bbn1(1,kkk)=ui(1,i+1);
%     
%     [m_bbn1,i_bbn1]=max(real(bbn1));
%    
%     bbn2(1,kkk)=ui(1,i+2);
%     
%     [m_bbn2,i_bbn2]=max(real(bbn2));
% 
%     bbn3(1,kkk)=ui(1,i+3);
%     
%     [m_bbn3,i_bbn3]=max(real(bbn3));
%     
%     bbn4(1,kkk)=ui(1,i+4);
%     
%     [m_bbn4,i_bbn4]=max(real(bbn4));
   



sim_phase=length(bbn);


for i12=1:length_1
arr5(i12)=phase1(actt1(i12));

end

% 
% y13(1)=0;
% y14(1)=0;
% 
% k5=2;
% for i=1:5
%     
%    for j=1:16
%        
%       y133(i,j)=y13(k5)+y13(k5-1);
%       y144(i,j)=y14(k5)+y14(k5-1);
%    end
%    
%    k5=k5+1;
% end
%    
%     
%     
%     
clear i6
clear j6
% 
% for x666=1:5
%     for j6=1:4
%          x7(j6,x666)=cell{j6}(x666);
%          x6(j6)=cell{j6}(x666);
%     end
%         l6=max(real(x666));
%         l7=find(x666==l6);
%         l1(x666)=l7;
%         x6=0;
% end
% for k11=1:5
%     f_ph(k11)=phase(l1(k11));
% end
for i9=1:length(arr5)
    states(i9)=atand((imag(arr5(i9)))/real(arr5(i9)));
    if(real(arr5(i9))<0 && imag(arr5(i9))>=0)
         states(i9)=(atand((imag(arr5(i9)))/(real(arr5(i9)))))+180;
    end
    if (states(i9)==-90)
        states(i9)=270;
    end
end


for i8=1:length(states)
    
    if (i8==1)
        
        if(states(i8)==90)
            rec_bits(i8)=1;
        else
            rec_bits(i8)=-1;
        end
    end
    
    
    
        if (i8<length_1+1 && i8~=1)
  if (states(i8)>states(i8-1))
      
     
      rec_bits(i8)=+1;
      
  else
      rec_bits(i8)=-1;
            
  end
        end
    
    if ( states(i8)==270 && i8~=1 && states(i8-1)==0 )
        if (states(i8)>states(i8-1))
        rec_bits(i8)=-1;
        else
            rec_bits(i8)=1;
        end
    end
    if ( states(i8)==0 && states(i8-1)==270 && i8~=1)
        if (states(i8)>states(i8-1))
        rec_bits(i8)=-1;
        else
            rec_bits(i8)=1;
        end
    end
end
% received=zeros(1,5);

% for j4=1:length(received)
%     
% received(j4)=(atand((imag(f_ph(j4)))/(real(f_ph(j4)))));
% 
%   else
%   end
% 
%   
% end
% for j=1:length(received)
%     if(real(f_ph(j))>=0 && imag(f_ph(j))>=0)
%     received(j)=(atand((imag(f_ph(j)))/(real(f_ph(j)))));
%     else if(real(f_ph(j))<=0 && imag(f_ph(j))>=0)
%             if((atand((imag(f_ph(j)))/(real(f_ph(j)))))>=90 && (atand((imag(f_ph(j)))/(real(f_ph(j)))))<=270)
%                 received(j)=(atand((imag(f_ph(j)))/(real(f_ph(j)))));
%             else
%                 a9=(atand(abs((imag(f_ph(j)))/(real(f_ph(j))))));
%                 b9=180-a9;
%                 received(j)=b9-360;
%             end
%     else if(real(f_ph(j))<=0 && imag(f_ph(j))<=0)
%         if((atand((imag(f_ph(j)))/(real(f_ph(j)))))>=90 && (atand((imag(f_ph(j)))/(real(f_ph(j)))))<=270)
%                 received(j)=(atand((imag(f_ph(j)))/(real(f_ph(j)))));
%             else
%                 a10=(atand(abs((imag(f_ph(j)))/(real(f_ph(j))))));
%                 b10=270-a10;
%                 received(j)=b10-360;
%             end
%     else if(real(f_ph(j))>=0 && imag(f_ph(j))<=0)
%            if((atand((imag(f_ph(j)))/(real(f_ph(j)))))>=90 && (atand((imag(f_ph(j)))/(real(f_ph(j)))))<=270)
%                 received(j)=(atand((imag(f_ph(j)))/(real(f_ph(j)))));
%             else
%                 a10=(atand(abs((imag(f_ph(j)))/(real(f_ph(j))))));
%                 b10=360-a10;
%                 received(j)=b10-360;
%             end
%         end
%         end
%         end
%     end
%end
subplot 311
stem(a1)
title ('INPUT BITS')
subplot 312
stem(rec_bits)
title ('RECEIVED BITS')
subplot 313
stem(transmitted1)
title ('CONTINOUS PHASE')

