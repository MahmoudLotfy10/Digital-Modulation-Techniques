clc;
clear;
close all;
%%
n_bits = (10^5)*24; %24=2*3*4 QPSK=>2 ,, 8PSK=>3 ,, 16QAM=>4
data = randi([0, 1], n_bits, 1);
Eb = 1;
snr_db = (-4:1:14);
img = 1i;

%%

%%%%%%%%%%%%%%%%%%%%%BPSK%%%%%%%%%%%%%%%%%

%%The Mapper
BPSK_Mapped =zeros([n_bits 1]);
for i=1:n_bits
    if(data(i,1)==1)
        BPSK_Mapped(i,1)=sqrt(Eb);
    else
         BPSK_Mapped(i,1)=-sqrt(Eb);
    end
end
%%The channel
BER_theoretical_BPSK = zeros(size(snr_db));
BER_BPSK = zeros(size(snr_db));
for j = 1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    Eb_BPSK = 1;
    No = Eb_BPSK / snr_linear;
    Eavg=mean(abs(BPSK_Mapped).^2);
    var = sqrt(sqrt(Eavg/1)*(No/2));
    AWGN_Channel_Noise_BPSK = var .* randn(n_bits, 1);
    BPSK_noise = BPSK_Mapped + AWGN_Channel_Noise_BPSK;
    
    %%Demapper
    BPSK_demapped = zeros(size(BPSK_noise));
    for i = 1:n_bits
        if BPSK_noise(i) > 0
            BPSK_demapped(i) = 1;
        else
            BPSK_demapped(i) = 0;
        end
    end
    BPSK_demapped=reshape(BPSK_demapped',[1 n_bits]);
    BPSK_demapped=BPSK_demapped';
    
%%BER
    N_error_BPSK=0;
    for i = 1:length(data)
        if((data(i) ==  1) && (BPSK_demapped(i) == 0))
            N_error_BPSK = N_error_BPSK+1;
        end
        if((data(i) == 0) && (BPSK_demapped(i) ==1 ))
            N_error_BPSK = N_error_BPSK+1;
        end
    end
    BER_theoretical_BPSK(j) = 0.5 * erfc(sqrt(Eb / No));
    BER_BPSK(j) = N_error_BPSK / n_bits;
end
figure;
plot(real(BPSK_Mapped), imag(BPSK_Mapped),'bo','linewidth',2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-4 4], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-2 2 -2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('BPSK Constellation ');
grid on;
axis equal;

figure('Name', 'BER of BPSK');
semilogy(snr_db, BER_theoretical_BPSK, 'b --');
hold on;
semilogy(snr_db, BER_BPSK, 'r -o');
hold off;
title('BER of BPSK');
xlim([-4  14]);
xlabel('Eb/No (dB)');
ylabel('BER');
ylim([10^-4  10]);
grid on;
legend('BPSK theoretical', 'BPSK Practical');
%%
%%%%%%%%%%%%%%%%%%%%%QPSK Gray Coding%%%%%%%%%%%%%%%%%
%%the Mapper
data_QPSK_Gray=reshape(data,[2,n_bits/2]);
data_QPSK_Gray=data_QPSK_Gray';
QPSK_Gray_Mapped=zeros([n_bits/2 1]);

for i=1:n_bits/2
    index=data_QPSK_Gray(i,:);
    if index==[1 0]
        QPSK_Gray_Mapped(i)=1-img;
    elseif index==[1 1]
        QPSK_Gray_Mapped(i)=1+img;
    elseif  index==[0 0]
        QPSK_Gray_Mapped(i)=-1-img;
    elseif  index==[0 1]
        QPSK_Gray_Mapped(i)=-1+img;
    end
end

BER_QPSK_Gray=zeros(size(snr_db));
BER_theoretical_QPSK_Gray=zeros(size(snr_db));
%%The channel stage
for j=1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    EB_QPSK = 1;
    No = EB_QPSK/snr_linear; 
    Eavg=mean(abs(QPSK_Gray_Mapped).^2);
    var = sqrt((Eavg/2)*(No/2));
    real_noise=var.* randn(1 , size(QPSK_Gray_Mapped,1));
    img_noise = var .* randn(1 , size(QPSK_Gray_Mapped,1));
    noise=(real_noise+(img*img_noise))';
    noise_QPSK_Gray = QPSK_Gray_Mapped + noise;
    
%%The Demapper stage    
QPSK_Gray_demapped=zeros([n_bits/2 2]);
    for i=1:n_bits/2
        real_number=real(noise_QPSK_Gray(i));
        imag_number=imag(QPSK_Gray_Mapped(i));
        if (real_number>=0) && (imag_number>=0)
            QPSK_Gray_demapped(i,:)=[1 1];
        elseif (real_number>=0) && (imag_number<=0)
            QPSK_Gray_demapped(i,:)=[1 0];
        elseif (real_number<=0) && (imag_number>=0)
            QPSK_Gray_demapped(i,:)=[0 1];
        elseif (real_number<=0) && (imag_number<=0)
            QPSK_Gray_demapped(i,:)=[0 0];
        end
    end
    QPSK_Gray_demapped=reshape(QPSK_Gray_demapped',[1 n_bits]);
    QPSK_Gray_demapped=QPSK_Gray_demapped';
    
%%BER 
N_bit_error_QPSK_gray=0;
   for i = 1:length(data)
        if((data(i) ==  1) && (QPSK_Gray_demapped(i) == 0))
            N_bit_error_QPSK_gray = N_bit_error_QPSK_gray+1;
        end
        if((data(i) == 0) && (QPSK_Gray_demapped(i) ==1 ))
            N_bit_error_QPSK_gray = N_bit_error_QPSK_gray+1;
        end
    end
    
    BER_QPSK_Gray(j)= 2*N_bit_error_QPSK_gray/length(data);
    BER_theoretical_QPSK_Gray(j)=0.5* erfc(sqrt(EB_QPSK/No));
end
figure;
plot(real(QPSK_Gray_Mapped), imag(QPSK_Gray_Mapped),'bo','linewidth',2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-4 4], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-2 2 -2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('QPSK Gray Constellation ');
grid on;
axis equal;

figure('Name' , 'BER of QPSK Gray');
semilogy(snr_db , BER_theoretical_QPSK_Gray,'b --');
hold on;
semilogy(snr_db , BER_QPSK_Gray,'r -o');
hold off; 
title('BER of QPSK Gray');
xlabel('EB/No(dB)');  
xlim([-4  14]);
ylabel('BER'); 
ylim([10^-4  10]);
grid on;
legend 'QPSK Gray theoretical' 'QPSK Gray Practical';
%% 
%%%%%%%%%%%%%%%%%%% QPSK Not Gray %%%%%%%%%%%%%%%
%%The Mapper
data_QPSK_notGray=reshape(data,[2,n_bits/2]);
data_QPSK_notGray=data_QPSK_notGray';
QPSK_Not_Gray_Mapped=zeros([n_bits/2 1]);

BER_QPSK_notGray=zeros(size(snr_db));
BER_theoretical_QPSK_notGray=zeros(size(snr_db));
for i=1:n_bits/2
    index=data_QPSK_notGray(i,:);
    if index==[1 0]
        QPSK_Not_Gray_Mapped(i)=1+img;
    elseif index==[1 1]
        QPSK_Not_Gray_Mapped(i)=1-img;
    elseif  index==[0 0]
        QPSK_Not_Gray_Mapped(i)=-1-img;
    elseif  index==[0 1]
        QPSK_Not_Gray_Mapped(i)=-1+img;
    end
end

%%The channel stage
for j=1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    EB_QPSK = 1;
    No = EB_QPSK/snr_linear; 
    Eavg=mean(abs(QPSK_Not_Gray_Mapped).^2);
    var = sqrt((Eavg/2)*(No/2)); 
    real_noise=var.* randn(1 , size(QPSK_Not_Gray_Mapped,1));
    img_noise = var .* randn(1 , size(QPSK_Not_Gray_Mapped,1));
    noise=(real_noise+(img*img_noise))';
    noise_QPSK_not_gray = QPSK_Not_Gray_Mapped +noise;
    
%%The Demapper stage
    demapped_QPSK_not_gary=zeros([n_bits/2 2]);
    for i=1:n_bits/2
        real_number=real(noise_QPSK_not_gray(i));
        imag_number=imag(noise_QPSK_not_gray(i));
        if (real_number>=0) && (imag_number>=0)
            demapped_QPSK_not_gary(i,:)=[1 0];
        elseif (real_number>=0) && (imag_number<=0)
            demapped_QPSK_not_gary(i,:)=[1 1];
        elseif (real_number<=0) && (imag_number>=0)
            demapped_QPSK_not_gary(i,:)=[0 1];
        elseif (real_number<=0) && (imag_number<=0)
            demapped_QPSK_not_gary(i,:)=[0 0];
        end
    end
    demapped_QPSK_not_gary=reshape(demapped_QPSK_not_gary',[1 n_bits]);
    demapped_QPSK_not_gary=demapped_QPSK_not_gary';
    
%%BER
    N_bit_error_QPSK_not_gray=0;
    for i = 1:length(data)
        if((data(i) ==  1) && (demapped_QPSK_not_gary(i) == 0))
            N_bit_error_QPSK_not_gray = N_bit_error_QPSK_not_gray+1;
        end
        if((data(i) == 0) && (demapped_QPSK_not_gary(i) ==1 ))
            N_bit_error_QPSK_not_gray = N_bit_error_QPSK_not_gray+1;
        end
    end
    BER_QPSK_notGray(j)= N_bit_error_QPSK_not_gray/length(data);
    BER_theoretical_QPSK_notGray(j)= 0.5*erfc(sqrt(EB_QPSK/No));
end
figure;
plot(real(QPSK_Not_Gray_Mapped), imag(QPSK_Not_Gray_Mapped),'bo','linewidth',2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-4 4], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-2 2 -2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('QPSK Not Gray Constellation ');
grid on;
axis equal;

figure('Name', 'BER of QPSK Not Gray');
semilogy(snr_db, BER_theoretical_QPSK_notGray, 'b');
hold on;
semilogy(snr_db, BER_QPSK_notGray, 'r -o');
hold off;
title('BER of QPSK Not Gray');
xlabel('EB/No(dB)');  
xlim([-4  14]);
ylabel('BER'); 
ylim([10^-4  10]);
grid on;
legend('QPSK Not Gray theoretical', 'QPSK Not Gray Practical');



figure('Name', 'BER of QPSK');
semilogy(snr_db, BER_QPSK_Gray, 'g -*');
hold on;
semilogy(snr_db, BER_QPSK_notGray, 'r -o');
hold off;
title('Comparison between Practical QPSK Gray and Not Gray');
xlabel('EB/No(dB)');  
xlim([-4  14]);
ylabel('BER'); 
ylim([10^-4  10]);
grid on;
legend('QPSK  Gray Practical', 'QPSK Not Gray Practical');

%% 
%%%%%%%%%%%%%%%%%%%%%  8PSK  %%%%%%%%%%%%%%%%%%%%
EB_M_8PSK=1;
data_8PSK=reshape(data,[3,n_bits/3]);
data_8PSK=data_8PSK';
M_8PSK_Mapped=zeros([n_bits/3 1]);

BER_8PSK =zeros(size(snr_db));
BER_theoretical_8PSK =zeros(size(snr_db));
%%the Mapper
for i=1:n_bits/3
    index= data_8PSK(i,:);
    if index==[0 0 0]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*0);
    elseif index==[0 0 1]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*1);
    elseif index==[0 1 1]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*2);
    elseif index==[0 1 0]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*3);
    elseif index==[1 1 0]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*4);
    elseif index==[1 1 1]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*5);
    elseif index==[1 0 1]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*6);
    elseif index==[1 0 0]
        M_8PSK_Mapped(i)=  sqrt(EB_M_8PSK)*exp(1i*(pi/4)*7);
    end
end

%%The channel stage
for j=1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    EB_M_8PSK = 1;
    No = EB_M_8PSK/(snr_linear); 
    Eavg=mean(abs(M_8PSK_Mapped).^2);
    var = sqrt((Eavg/3)*(No/2)); 
    real_noise=var.* randn(1 , size(M_8PSK_Mapped,1));
    img_noise = var .* randn(1 , size(M_8PSK_Mapped,1));
    noise=(real_noise+(img*img_noise))';
    M_8PSK_noise = M_8PSK_Mapped +noise;
    
%%The Demapper stage
    demapped_8PSK=zeros([n_bits/3 3]);
    for i=1:n_bits/3
        Angle=M_8PSK_noise(i);
        %%condition after || to sure correct demapped in case of the sympole shifted by 2pi
        
        if   (angle(Angle)>=-pi/8) && (angle(Angle)<=pi/8)||((angle(Angle)+2*pi>=-pi/8)&& (angle(Angle)+2*pi<=pi/8))
            demapped_8PSK(i,:)=[0 0 0];
        elseif  (angle(Angle)>=pi/8) && (angle(Angle)<=3*pi/8)||((angle(Angle)+2*pi>=pi/8)&& (angle(Angle)+2*pi<=3*pi/8))
            demapped_8PSK(i,:)=[0 0 1];
        elseif  (angle(Angle)>=3*pi/8) && (angle(Angle)<=5*pi/8)||((angle(Angle)+2*pi>=3*pi/8)&& (angle(Angle)+2*pi<=5*pi/8))
            demapped_8PSK(i,:)=[0 1 1];
        elseif  (angle(Angle)>=5*pi/8) && (angle(Angle)<=7*pi/8)||((angle(Angle)+2*pi>=5*pi/8)&& (angle(Angle)+2*pi<=7*pi/8))
            demapped_8PSK(i,:)=[0 1 0];
        elseif  ((angle(Angle)>=7*pi/8) && (angle(Angle)<=9*pi/8))||((angle(Angle)+2*pi>=7*pi/8)&& (angle(Angle)+2*pi<=9*pi/8))
            demapped_8PSK(i,:)=[1 1 0];
        elseif  (angle(Angle)+2*pi>=9*pi/8) && (angle(Angle)+2*pi<=11*pi/8)||((angle(Angle)+2*pi>=9*pi/8)&& (angle(Angle)+2*pi<=11*pi/8))
            demapped_8PSK(i,:)=[1 1 1];
        elseif  (angle(Angle)+2*pi>=11*pi/8) && (angle(Angle)+2*pi<=13*pi/8)||((angle(Angle)+2*pi>=11*pi/8)&& (angle(Angle)+2*pi<=13*pi/8))
            demapped_8PSK(i,:)=[1 0 1];
        elseif  (angle(Angle)+2*pi>=13*pi/8) && (angle(Angle)<=15*pi/8)||((angle(Angle)+2*pi>=13*pi/8)&& (angle(Angle)+2*pi<=15*pi/8))
            demapped_8PSK(i,:)=[1 0 0];
        end
    end
    demapped_8PSK=reshape(demapped_8PSK',[1 n_bits]);
    demapped_8PSK=demapped_8PSK';

%%BER
    N_bit_error_8PSK=0;
    for i=1:length(data)
        if((data(i) ==  1) && (demapped_8PSK(i) == 0))
            N_bit_error_8PSK = N_bit_error_8PSK+1;
        end
        if((data(i) == 0) && (demapped_8PSK(i) ==1 ))
            N_bit_error_8PSK = N_bit_error_8PSK+1;
        end
    end
    BER_8PSK(j) =N_bit_error_8PSK/length(data);
    BER_theoretical_8PSK(j) = erfc(sqrt(3/(No))*sin(pi/8))/3;
end
figure;
plot(real(M_8PSK_Mapped), imag(M_8PSK_Mapped),'bo','linewidth',2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-4 4], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-2 2 -2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('8PSK Constellation ');
grid on;
axis equal;


figure('Name' , 'BER of 8PSK');
semilogy( snr_db,BER_theoretical_8PSK,'b--');
hold on;
semilogy( snr_db,BER_8PSK,'r -o');
hold off;
title('BER of 8PSK');
xlabel('EB/No(dB)');  
xlim([-4  14]);
ylabel('BER'); 
ylim([10^-4  10]);
grid on;
legend('8PSK theoretical','8PSK Practical')
%%
%%%%%%%%%%%%%%%%%%% 16QAM %%%%%%%%%%%%%%%
%%The Mapper
EB_16_QAM=1;
data_16QAM=reshape(data,[4,n_bits/4]);
data_16QAM=data_16QAM';
QAM_16_Mapped=zeros([n_bits/4 1]);

BER_16_QAM=zeros(size(snr_db));
BER_theoretical_16_QAM=zeros(size(snr_db));
for i=1:n_bits/4
    index=data_16QAM(i,:);
    if index==[1 1 1 1]
        QAM_16_Mapped(i)=1+img;
    elseif index==[1 0 1 1]
        QAM_16_Mapped(i)=3+img;
    elseif  index==[1 1 1 0]
        QAM_16_Mapped(i)=1+3*img;
    elseif  index==[1 0 1 0]
        QAM_16_Mapped(i)=3+3*img;
        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif  index==[1 1 0 1]
        QAM_16_Mapped(i)=1-img;
    elseif  index==[1 0 0 1]
        QAM_16_Mapped(i)=3-img;
    elseif  index==[1 1 0 0]
        QAM_16_Mapped(i)=1-3*img;
    elseif  index==[1 0 0 0]
        QAM_16_Mapped(i)=3-3*img;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     elseif  index==[0 1 1 1]
        QAM_16_Mapped(i)=-1+img;
     elseif  index==[0 1 1 0]
        QAM_16_Mapped(i)=-1+3*img;
     elseif  index==[0 0 1 1]
        QAM_16_Mapped(i)=-3+1*img;
     elseif  index==[0 0 1 0]
        QAM_16_Mapped(i)=-3+3*img;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     elseif  index==[0 1 0 1]
        QAM_16_Mapped(i)=-1-1*img;
     elseif  index==[0 0 0 1]
        QAM_16_Mapped(i)=-3-1*img;
     elseif  index==[0 1 0 0]
        QAM_16_Mapped(i)=-1-3*img;
     elseif  index==[0 0 0 0]
        QAM_16_Mapped(i)=-3-3*img;
    end
end

for j=1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    EB_16QAM = 1;
    No = EB_16QAM/snr_linear; 
    Eavg=mean(abs(QAM_16_Mapped).^2);
    var = sqrt((Eavg/4)*(No/2)); 
    real_noise=var.* randn(1 , size(QAM_16_Mapped,1));
    img_noise = var .* randn(1 , size(QAM_16_Mapped,1));
    noise=(real_noise+(img*img_noise))';
    noise_16_QAM = QAM_16_Mapped +noise;
    
%%The Demapper stage
   demapped_16_QAM=zeros([n_bits/4 4]);
    for i=1:n_bits/4
        real_number=real(noise_16_QAM(i));
        imag_number=imag(noise_16_QAM(i));
        
        
        if (real_number>=0)&& (real_number<=2)
            if(imag_number>=0)&&(imag_number<=2)
                demapped_16_QAM(i,:)=[1 1 1 1];
            elseif (imag_number<=0)&&(imag_number>=-2)
                demapped_16_QAM(i,:)=[1 1 0 1];
            elseif (imag_number>=2)
                demapped_16_QAM(i,:)=[1 1 1 0];
            elseif (imag_number<=-2)
                demapped_16_QAM(i,:)=[1 1 0 0];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif (real_number>=2)
            if(imag_number>=0)&&(imag_number<=2)
                demapped_16_QAM(i,:)=[1 0 1 1];
            elseif (imag_number<=0)&&(imag_number>=-2)
                demapped_16_QAM(i,:)=[1 0 0 1];
            elseif (imag_number>=2)
                demapped_16_QAM(i,:)=[1 0 1 0];
            elseif (imag_number<=-2)
                demapped_16_QAM(i,:)=[1 0 0 0];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            elseif (real_number<=0)&&(real_number>=-2)
            if(imag_number>=0)&&(imag_number<=2)
                demapped_16_QAM(i,:)=[0 1 1 1];
            elseif (imag_number<=0)&&(imag_number>=-2)
                demapped_16_QAM(i,:)=[0 1 0 1];
            elseif (imag_number>=2)
                demapped_16_QAM(i,:)=[0 1 1 0];
            elseif (imag_number<=-2)
                demapped_16_QAM(i,:)=[0 1 0 0];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (real_number<=-2)
            if(imag_number>=0)&&(imag_number<=2)
                demapped_16_QAM(i,:)=[0 0 1 1];
            elseif (imag_number<=0)&&(imag_number>=-2)
                demapped_16_QAM(i,:)=[0 0 0 1];
            elseif (imag_number>=2)
                demapped_16_QAM(i,:)=[0 0 1 0];
            elseif (imag_number<=-2)
                demapped_16_QAM(i,:)=[0 0 0 0];
            end
            
            
        end
    end
    demapped_16_QAM=reshape(demapped_16_QAM',[1 n_bits]);
    demapped_16_QAM=demapped_16_QAM';
    
%%BER
    N_bit_error_16_QAM=0;
    for i = 1:length(data)
        if((data(i) ==  1) && (demapped_16_QAM(i) == 0))
            N_bit_error_16_QAM = N_bit_error_16_QAM+1;
        end
        if((data(i) == 0) && (demapped_16_QAM(i) ==1 ))
            N_bit_error_16_QAM = N_bit_error_16_QAM+1;
        end
    end
    BER_16_QAM(j)= N_bit_error_16_QAM/length(data);
    BER_theoretical_16_QAM(j)= (1.5/4)*erfc(sqrt(1/(2.5*No)));
end
figure;
plot(real(QAM_16_Mapped), imag(QAM_16_Mapped),'bo','LineWidth', 2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-6 6], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-4 4 -4 4]);
xlabel('In-Phase');
ylabel('Quadrature');
title('16QAM Constellation ');
grid on;
axis equal;

figure('Name' , 'BER of 16QAM');
semilogy( snr_db,BER_theoretical_16_QAM,'b--');
hold on;
semilogy( snr_db,BER_16_QAM,'r -o');
hold off;
title('BER of 16QAM');
xlabel('EB/No(dB)');  
xlim([-4  14]);
ylabel('BER'); 
ylim([10^-4  10]);
grid on;
legend('16QAM theoretical','16QAM Practical')
%%
%%%%%%%%%%%%%%%%%%%%%BFSK%%%%%%%%%%%%%%%%%
%%The Mapper Stage
EB_BFSK=1;
BFSK_Mapped=zeros([n_bits 1]);

BER_BFSK=zeros(size(snr_db));
BER_theoretical_BFSK=zeros(size(snr_db));

for i=1:n_bits
    if(data(i,1)==0)
        BFSK_Mapped(i,1)=sqrt(EB_BFSK);
    else
         BFSK_Mapped(i,1)=sqrt(EB_BFSK)*img;
    end
end
%%the channel stage
for j=1:length(snr_db)
    snr_linear = 10^(snr_db(j)/10);
    No = EB_BFSK/snr_linear; 
    Eavg=mean(abs(BFSK_Mapped).^2);
    var = sqrt((Eavg/1)*(No/2)); 
    real_noise=var.* randn(1 , size(BFSK_Mapped,1));
    img_noise = var .* randn(1 , size(BFSK_Mapped,1));
    noise=(real_noise+(img*img_noise))';
    BFSK_noise = BFSK_Mapped +noise;
   
    %%the demapper
    
   BFSK_demapped=zeros([n_bits 1]);
   for i=1:n_bits
        Angle=angle(BFSK_noise(i));
        if (Angle<=pi/4)||((Angle+(2*pi))<=pi/4)
            BFSK_demapped(i,:)=0;
        else
            BFSK_demapped(i,:)=1;
        end
   end
    BFSK_demapped=reshape(BFSK_demapped',[1 n_bits]);
    BFSK_demapped=BFSK_demapped';
    
    %%BER
    N_bit_error_BFSK=0;
    for i=1:length(data)
        if((data(i) ==  1) && (BFSK_demapped(i) == 0))
            N_bit_error_BFSK = N_bit_error_BFSK+1;
        end
        if((data(i) == 0) && (BFSK_demapped(i) ==1 ))
            N_bit_error_BFSK = N_bit_error_BFSK+1;
        end
    end
    BER_BFSK(j) =N_bit_error_BFSK/length(data);
    BER_theoretical_BFSK(j) = 0.5*erfc(sqrt(Eb / (2*No)));
end
figure;
plot(real(BFSK_Mapped), imag(BFSK_Mapped),'bo','LineWidth', 2);
hold on;
line([0 0], [-4 4], 'Color', 'k', 'LineStyle', '-'); 
line([-4 4], [0 0], 'Color', 'k', 'LineStyle', '-'); 
hold off;
axis([-2 2 -2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('BFSK Constellation ');
grid on;
axis equal;

figure('Name', 'BER of BFSK');
semilogy(snr_db, BER_theoretical_BFSK, 'b --');
hold on;
semilogy(snr_db, BER_BFSK, 'r -o');
hold off;
title('BER of BFSK');
xlim([-4  14]);
xlabel('Eb/No (dB)');
ylabel('BER');
ylim([10^-4  10]);
grid on;
legend('BFSK theoretical', 'BFSK Practical');


%%%%%%%plot all Theoritical%%%%%%%
figure('Name', 'All Theoritical BER');
semilogy(snr_db, BER_theoretical_BPSK);
hold on;
semilogy(snr_db , BER_theoretical_QPSK_Gray);
semilogy(snr_db, BER_theoretical_QPSK_notGray);
semilogy( snr_db,BER_theoretical_8PSK);
semilogy( snr_db,BER_theoretical_16_QAM);
semilogy(snr_db, BER_theoretical_BFSK);
hold off;
title('All Theoritical BER');
xlim([-4  14]);
xlabel('Eb/No (dB)');
ylabel('BER');
ylim([10^-4  10]);
grid on;
legend('BPSK theoretical', 'QPSK Gray theoretical','QPSK Not Gray theoretical','BER 8PSK theoretical ','BER 16QAM theoretical','BER BFSK theoretical');
%%%%%%%plot all Practical%%%%%%%
figure('Name', 'All Practical BER');
semilogy(snr_db, BER_BPSK);
hold on;
semilogy(snr_db , BER_QPSK_Gray);
semilogy(snr_db, BER_QPSK_notGray);
semilogy( snr_db,BER_8PSK);
semilogy( snr_db,BER_16_QAM);
semilogy(snr_db, BER_BFSK);
hold off;
title('All Practical BER');
xlim([-4  14]);
xlabel('Eb/No (dB)');
ylabel('BER');
ylim([10^-4  10]);
grid on;
legend('BPSK ', 'QPSK Gray ','QPSK Not Gray ','BER 8PSK  ','BER 16QAM ','BER BFSK ');

%%%%%%%plot all %%%%%%%
figure('Name', 'All Practical and Theoritical BER');
semilogy(snr_db, BER_theoretical_BPSK);
hold on;
semilogy(snr_db, BER_BPSK);
semilogy(snr_db , BER_theoretical_QPSK_Gray);
semilogy(snr_db , BER_QPSK_Gray);
semilogy(snr_db, BER_theoretical_QPSK_notGray);
semilogy(snr_db, BER_QPSK_notGray);
semilogy( snr_db,BER_theoretical_8PSK);
semilogy( snr_db,BER_8PSK);
semilogy( snr_db,BER_theoretical_16_QAM);
semilogy( snr_db,BER_16_QAM);
semilogy(snr_db, BER_theoretical_BFSK);
semilogy(snr_db, BER_BFSK);
hold off;
title('All Practical and Theoritical BER');
xlim([-4  14]);
xlabel('Eb/No (dB)');
ylabel('BER');
ylim([10^-4  10]);
grid on;
legend('BPSK theoretical','BPSK ',  'QPSK Gray theoretical','QPSK Gray ','QPSK Not Gray theoretical','QPSK Not Gray ','BER 8PSK theoretical ','BER BFSK theoretical','BER 8PSK  ','BER 16QAM theoretical','BER 16QAM ','BER BFSK ');


%%%PSD
Eb = 1;
Tb = 0.07;
delta_f = 1/Tb;
f1 = 1/Tb;
f2 = f1 + delta_f;
fs = 100;
t_tot = 7;                      
t = 0:1/fs:Tb; 
samples = 7;
num_bits = 100;
ensemble_size = 500;
realizations=samples*num_bits;
f = (-realizations/2 : realizations/2 - 1) * fs / realizations; 
x0 = 0.5 / Tb;
l=hamming(realizations+7);
data = randi([0,1], ensemble_size, num_bits + 1 );
data = repelem(data, 1, 7);
Mapped_data = zeros(ensemble_size,(realizations+7));
c=0;
for i = 1:ensemble_size
    for j = 1:realizations
        if(data(i,j)==1)
            Mapped_data(i,j) = sqrt(2*Eb/Tb);  %S1BB
        else
            if c<7
                c=c+1;
            Mapped_data(i,j)= sqrt(2*Eb/Tb)*(cos(2*pi*delta_f*t(1,c))+1i*sin(2*pi*delta_f*t(1,c))); %S2BB
            end
            if c==7
                c=0;
            end
        end
    end
end
Mapped_data = Mapped_data.*l';
delay = randi ([0,samples-1], ensemble_size, 1 );
Mapped_data_delayed = zeros(realizations,(samples*num_bits)); 

for i = 1:ensemble_size
    Mapped_data_delayed(i, :) = Mapped_data(i, 1 + delay(i) : realizations + delay(i));
end
ACF = zeros (1,(realizations)) ;


    for tau = -349:350
        
        prodofconj = conj(Mapped_data_delayed(:, realizations/2)) .* Mapped_data_delayed(:, realizations/2 + tau);
       
        ACF(tau + realizations/2) = sum(prodofconj) / realizations;
    end

PSD  =abs(fftshift(fft(ACF)));
PSD_normalized = (PSD)/fs ;
figure;
plot(f*Tb, PSD_normalized);
xlabel('Freq (HZ)');
ylabel('Amplitude');
xlim([-1 2]);
ylim([0 2]);
grid on;
title('Practical Psd');

delta1 = 99 * double(abs(f - (x0)) < 0.01); 
delta2 =99 * double(abs(f - (-x0)) < 0.01); 
theoretical_psd = ((2/Tb) * (delta1 + delta2)) +((8 * cos(pi * Tb * f).^2) ./ (pi^2 * (4 * Tb^2 * f.^2 - 1).^2));
figure;
plot((f*Tb)+.5, (theoretical_psd));
xlabel('Frequency');
ylabel('Amplitude');
title('Theoretical psd');
grid on;
ylim([0 2]);
xlim([-2 2]);

figure;
plot((f*Tb)+.5, theoretical_psd);
hold on;
plot(f*Tb, PSD_normalized,'r');
hold off;
xlabel('Freq (HZ)');
ylabel('Amplitude');
title('Comparision between Theoretical and Practical psd');
grid on;
legend 'theoretical psd' 'Practical psd';
ylim([0 2]);
xlim([-2 2]);
