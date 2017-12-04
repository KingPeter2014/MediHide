% Copyright, Eze P.U., Udaya P. 2017.
% Email: peze@student.unimelb.edu.au
% Algorithm designed by Eze P.U
clear;
clc;
blocksize=8;
desiredcorrelation= 2;
compressionRatio = 2;

goldcode=[-1, 1, 1, -1, 1, -1, 1, 1;
    1, -1, -1, 1, 1, 1, 1, -1;
    1, 1, 1, 1, 1, -1, -1, -1;
    -1, -1, -1, 1, -1, 1, -1, 1;
    -1, -1, 1, 1, -1, -1, 1, -1;
    -1, -1, 1, -1, -1, 1, -1, 1;
    1, -1, 1, 1, -1, -1, -1, 1;
    1, 1, -1, 1, -1, -1, -1,1];



dicomlist = dir(fullfile(pwd,'DataSets/','*.dcm'));
for cnt = 1 : numel(dicomlist)
    I{cnt} = dicomread(fullfile(pwd,'DataSets/',dicomlist(cnt).name));  
end

% =======DEFINE GLOBAL VARIABLES FO RUNNING THE SIMULATION=========
withoutwatermark=[];
zeroAlpha = [];
oneAlpha=[];
embeddable =0;alphamax = 25; dis = [];
extractzerocorrelation=[];% correlations to extract watermark in blocks
onesextracted=[];extractedwatermark='';
ssmglobal=[];psnglobal=[];
ssm0=[];psn0=[];ssm1=[];psn1=[];
BER=0;berOverall=[];fnegative=[];bervalue = [];

entropyglobal=[];rentropyglobal=[];
entropylocal=[];rentropylocal=[];

entropyglobalstego = [];entropylocalstego = [];

overallTotalextracted=0;grandTotalEmbed=0;
epsilon = 0.5;% Tolerance for overflow, underflow and unintentional attack

for cnt = 1 : numel(dicomlist)
    dicom=strcat('DataSets/',dicomlist(cnt).name);
    % cell2mat()
     info=dicominfo(dicom);% extract the metadata from the image
     metadata = info;
     im=dicomread(info);
    %  im = I(1);
     pixeldepth=info.BitDepth;
    peakval=2^(pixeldepth);
     if length(size(im)) >2
        im=rgb2gray(im);
        display('Coverted to gray successfully')
     end
    [row,col]=size(im);
   
    im_uncast = im;%NOT cast from int16 to uint16
     im = uint16(im); % Necessary to ensure pixel values are not accepted at underflow and overflow values
    
    block=im(1:blocksize,1:blocksize);
    [N, M]=size(block);
    message=round(rand(row/blocksize,col/blocksize));% Experimental message from random bits
%     load message.mat
    %===============================C4S ALGORITHMIC ===============================
    stringWatermark = binaryToStringBits(message);
    length_of_watermark = length(stringWatermark);
    fractional = mod(length_of_watermark,compressionRatio);
    if fractional ~=0 
       padding =  compressionRatio - fractional;
       stringWatermark = strcat(stringWatermark,dec2bin(0,padding)); % pad with zeros
    end
    embedding_Channels = 2^compressionRatio;
    channel_gap = 4*epsilon;
    nearNegZeroCorrBitGroup = (embedding_Channels/2)-1; % Corresponds to correlation of -channel_gap
    nearPosZeroCorrBitGroup = (embedding_Channels/2);% Corresponds to correlation of channel_gap
    startindex = 1;
    %========================END C4S================================
    i=1;j=1;
    caption=strcat('Watermark Embedding and Extraction for:', dicom,'. Please wait...');
    display(caption);
    bitToEmbed=1;
    thisextractedwatermark=[];% Only for current image
    finalstego = zeros(row,col);
    finalstego=im;
    withTampterDetected = zeros(row,col);
    modifiedcover = zeros(row,col);
    totalmodifiedpixels=0;
    none=0;
    x=1;y=1;
    
    ber=0; % error bits
    falsenegative =0;
    while i<row && startindex <= length_of_watermark
        while j<col && startindex <= length_of_watermark
            block=im(i:i+7,j:j+7);
            withoutwatermark=[withoutwatermark linearCorrelation(block,goldcode)];
            
            %===============================C4S ALGORITHM ===============================
            thisWatermarkGroup = substr(stringWatermark,startindex,compressionRatio);
            watermarkGroupBitValue = bin2dec(thisWatermarkGroup);
            caption=strcat('Embedded :', thisWatermarkGroup);
%             display(caption)
            signeddesiredcorrelation = computeDesiredCorrelation( watermarkGroupBitValue,compressionRatio,epsilon);
            desiredcorrelation = abs(signeddesiredcorrelation);
%             return;
            %========================END C4S================================
            if sign(signeddesiredcorrelation) == -1
                bitToEmbed=1; % Negative correlation bitgroups
            else
                bitToEmbed=0;% Positive correlation bitgroups
            end 
                
            alpha=computeAlpha(block,desiredcorrelation,bitToEmbed,goldcode,'ASS');

            %==== HISTOGRAM EQUALIZATION TO OVERCOME OVERFLOW AND UNDERFLOW
            [modpixels,block]=HistogramEqualize(block,pixeldepth,2*alpha);%Equalise histogram to remove underflow and overflow
            modifiedcover(i:i+7,j:j+7)=block;
            totalmodifiedpixels=totalmodifiedpixels + modpixels;
            alpha=computeAlpha(block,desiredcorrelation,bitToEmbed,goldcode,'ASS');
            %================Ends Here===================%
            [ answer,hsi,di ] = isC4SEmbeddable(block,goldcode,alphamax, desiredcorrelation );
            if answer == 1
                dis = [dis di];
                embeddable = embeddable +1;
            end
            zeroAlpha=[zeroAlpha alpha];
            watermark=(alpha.*goldcode);
            if (bitToEmbed==0) % Embed positive correlation bitgroup watermark
                stego=(double(block)+ watermark);
            else
                stego=(double(block)- watermark);% Embed negative correlation bitgroup watermark
               
            end

            SSM = ssim(uint16(stego),block);% more than 8-bit depth
            PSN = psnr(uint16(stego),block,peakval);
            
            entropylocal = [entropylocal entropy(block)];
            entropylocalstego =[entropylocalstego entropy(uint16(stego))];
            retrylocal = relativeEntropy(uint16(stego),block);
    %         PSN = psnr(uint8(stego),block);
    %         SSM = ssim(uint8(stego),block);
            ssm0=[ssm0 SSM];psn0=[psn0 PSN];
            
            rentropylocal=[rentropylocal retrylocal];
  
    %============================ATTACKS=================================
%     blur_gaussian_filter = imgaussfilt(uint16(stego),2);% Gaussian Filter with sigma = 2 
%     guassian_noise=imnoise(uint16(stego),'gaussian',0,0.005); % Guassian noise, mean and var
%     spepper_noise = imnoise(uint16(stego),'salt & pepper',0.05);
%     speckle_noise=imnoise(uint16(stego),'speckle',0.05);
%     poisson_noise=imnoise(uint16(stego),'poisson');
%     cropped=imcrop(uint16(stego));% Use crop tool to crop image
%     rotated=imrotate(uint16(stego),20); % Rotate image by 20%
%       J = imresize(uint16(stego), 0.5);
%       K = medfilt2(uint16(stego); % Median Filter with 3x3
%       h = ones(3,3) / 9;
% %     avg_filtered = imfilter(uint16(stego),h); % Average Filter
% %      compressed = imwrite(uint16(stego),'compressed.jpg','jpg','Comment','Compression attack');
%     histequalised = histeq(uint16(stego)); % Histogram Equalisation
%     contrastadjusted = adapthisteq(uint16(stego));
%     imadjust(uint16(stego)) % Adjust image intensity values
%      
    

%     =============EXTRACTION BEGINS HERE==========================
            corr0=linearCorrelation(uint16(stego),goldcode);
    %         corr0=linearCorrelation(uint8(stego),goldcode);
            finalstego(i:i+7,j:j+7)=stego;
            withTampterDetected(i:i+7,j:j+7)=stego;
            extractzerocorrelation= [extractzerocorrelation corr0];
           
    

            % C4S grouped-bit extraction and tamper detection
            groupbit = C4SDetectWatermark(corr0,compressionRatio,epsilon);
            
            
            if strcmp('0',groupbit) && compressionRatio ~=1 % No valid watermark extracted in the block
                groupbit = dec2bin(0,compressionRatio);
                extractedwatermark= strcat(extractedwatermark,groupbit);%False negative value - non-binary
                falsenegative =falsenegative+compressionRatio;
                withTampterDetected(i:i+7,j:j+7)=2^(pixeldepth);% Flagging potentially tampered blocks
                thisextractedwatermark = strcat(thisextractedwatermark,groupbit);
                display('No valid watermark extracted from this block') 
            elseif strcmp('9',groupbit)&& compressionRatio ==1
                extractedwatermark= strcat(extractedwatermark,groupbit);%False negative value - non-binary
                falsenegative =falsenegative+compressionRatio;
                withTampterDetected(i:i+7,j:j+7)=2^(pixeldepth);% Flagging potentially tampered blocks
                thisextractedwatermark = strcat(thisextractedwatermark,groupbit);
                display('No valid watermark extracted from this block') 
                
            else 
                extractedwatermark= strcat(extractedwatermark,groupbit);
                thisextractedwatermark = strcat(thisextractedwatermark,groupbit);
            end
            
            caption=strcat('Extracted :', groupbit);
%             display(caption)

            startindex = startindex + compressionRatio; % Go to next message group
            j=j+8;%Go to next column
        end
        j=1;
        i=i+8;% Go to next row
    end
    grandTotalEmbed = grandTotalEmbed + length_of_watermark;
    
    SSM = ssim(uint16(stego),block);% more than 8-bit depth     
    ssmglobal=[ssmglobal ssim(uint16(finalstego),im)];
    
    entropyglobal = [entropyglobal entropy(im)];
    entropyglobalstego = [entropyglobalstego entropy(uint16(finalstego))];
    
    psnglobal=[psnglobal psnr(uint16(finalstego),im,peakval)];
    
    rentropyglobal= [rentropyglobal relativeEntropy(uint16(finalstego),im)];
    
    ber = computeBER(stringWatermark,thisextractedwatermark);
    bervalue=[bervalue ber];
    BER = (ber)/length_of_watermark;
    berOverall=[berOverall BER];fnegative=[fnegative falsenegative];
    ber=0;falsenegative=0;BER=0;
    
    
    % ===========SAVE WATERMARKED AND POTENTIALLY TAMPERED IMAGES======
%     filename=strcat('c4s\StegoData\9bit\watermarkedDicom',num2str(cnt),'.dcm');
%     dicomwrite(uint16(finalstego),filename, metadata);
%     filename=strcat('c4s\StegoData\9bit\Tampered',num2str(cnt),'.dcm');
%     dicomwrite(uint16(withTampterDetected),filename, metadata);

    falseNegative = (falsenegative)/length_of_watermark;
    %dicomwrite(X, 'ct_copy.dcm', metadata, 'CreateMode', 'copy');
     

end % End for EACH image

ber_average = sum(bervalue)/grandTotalEmbed
% bervalue
% berOverall
averageGlobalPSNR = mean(psnglobal)

% Overall unwatermarked statistics
 meanTg = mean(withoutwatermark)% Mean Correlation without embedded data
%
 variancethreshold=var(withoutwatermark);
 standarddev= std(withoutwatermark);
 recommendedConstantCorrelation=meanTg+standarddev; % use this as p



figure
plot(zeroAlpha,psn0,'*')
xlabel('Embedding Strength,a')
ylabel('PSNR (dB)')
caption=strcat('Embedding Strength, a Vs PSNR at Compression: x',num2str(compressionRatio));
title( caption)
grid on
grid minor

% figure
% plot(zeroAlpha,ssm0,'o')
% xlabel('Embedding Strength,a')
% ylabel('SSIM')
% caption=strcat('Embedding Strength, a Vs SSIM at Compression: x',num2str(compressionRatio));
% title(caption)
% grid on
% grid minor

figure
bins=100;
histogram(withoutwatermark,bins)
% axis([-bins/2 bins/2 0 1500])
title('Host Signal Interference (HSI) for Image blocks')
xlabel('Linear Correlation')
ylabel('Number of Image Blocks')
grid on
grid minor


figure
bins=1000;
histogram(extractzerocorrelation,bins)

caption=strcat('Deviations from constant correlation, Compression: x',num2str(compressionRatio));
title(caption)
xlabel('Linear Correlation')
ylabel('Number of Image Blocks')
grid on
grid minor

figure
bins=30;
histogram(zeroAlpha,bins)

caption=strcat('Distribution of Embedding Strength, a at Compression: x',num2str(compressionRatio));
title(caption)
xlabel('Embedding strength,a')
ylabel('Number of Image Blocks')
grid on
grid minor

% figure
% %bins=30;
% histogram(ssm0,bins)
% caption=strcat('Distribution of SSIM at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('SSIM')
% ylabel('Number of Image Blocks')
% grid on
% grid minor

% figure
% %bins=30;
% histogram(psn0,bins)
% caption=strcat('Distribution of PSNR at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('PSNR(dB)')
% ylabel('Number of Image Blocks')
% grid on
% grid minor

% figure
% histogram(berOverall)
% caption=strcat('Distribution of BER at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('BER')
% ylabel('Number of Images')
% grid on
% grid minor

% figure
% plot(1:1:cnt,berOverall,'-*','LineWidth',2);
% xlabel('Image ID');ylabel('BER'); 
% caption=strcat('Image ID vs BER values at Compression: x',num2str(compressionRatio));
% title(caption)
% grid on
% grid minor

% figure
% histogram(bervalue)
% caption=strcat('Distribution of BER at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('BER (No. of Bits)')
% ylabel('Number of Images')
% grid on
% grid minor

% figure
% bins=20;
% histogram(ssmglobal)
% caption=strcat('Distribution of Global SSIM at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('Global SSIM Values')
% ylabel('Number of Images')
% grid on
% grid minor

% figure
% histogram(psnglobal)
% caption=strcat('Distribution of Global PSNR at Compression: x',num2str(compressionRatio));
% title(caption)
% xlabel('Global PSNR Values')
% ylabel('Number of Images')
% grid on
% grid minor

% figure
% plot(1:1:cnt,berOverall,'-*','LineWidth',2);
% xlabel('Image ID');ylabel('BER'); 
% caption=strcat('Image ID vs BER values at Compression: x',num2str(compressionRatio));
% title(caption)
% grid on
% grid minor

% figure 
% avgPSNR = [72.3974 74.3764 74.9794 74.4336 72.2662 68.7166 63.9252 58.3196 52.2041 46.0337];
% plot(1:1:10,avgPSNR,'-o','LineWidth',2);
% xlabel('Compression Rate');ylabel('Average PSNR(dB)'); 
% caption=strcat('Compression Rate Vs Average PSNR');
% title(caption)
% grid on
% grid minor

% figure
% plot(1:1:cnt,psnglobal,'*','LineWidth',2);
% xlabel('Image ID');ylabel('PSNR'); 
% caption=strcat('Image ID vs Global PSNR at Compression: x',num2str(compressionRatio));
% title(caption)
% grid on
% grid minor

% figure
% plot(rentropyglobal,psnglobal,'*','LineWidth',2);
% xlabel('Relative Entropy');ylabel('PSNR'); 
% caption=strcat('Relative Entropy vs PSNR at Compression: x',num2str(compressionRatio));
% title(caption)
% grid on
% grid minor

figure
plot(1:1:cnt,entropyglobal,'-r',1:1:cnt,entropyglobalstego,'-b');
xlabel('Image ID');ylabel('Global Image Entropy'); 
legend('Original','Stego')
caption=strcat('Image ID vs Original and Stego Global Entropy: x',num2str(compressionRatio));
title(caption)
grid on
grid minor

diff = [];
for i=1:length(entropylocal)
   diff1= entropylocal(1,i) - entropylocalstego(1,i);
    diff = [diff  diff1];
end
diff = sort(abs(diff));

figure
plot(1:1:length(entropylocal),entropylocal,'ob',1:1:length(entropylocal),entropylocalstego,'or');
xlabel('Image ID');ylabel('Local Block Entropy'); 
legend('Original','Stego')
caption=strcat('Image ID vs Original and Stego Local Entropy: x',num2str(compressionRatio));
title(caption)
grid on
grid minor

figure
plot(1:1:length(diff),diff,'-b','LineWidth',2);
xlabel('Image Sub-Blocks');ylabel('Entropy Difference'); 
caption=strcat('Image Sub-blocks vs Entropy Difference: x',num2str(compressionRatio));
title(caption)
grid on
grid minor

%=======Compare Entropy Differences==========

load entropydiff2.mat;
diff2=diff;
load entropydiff3.mat;
diff3=diff;
load entropydiff4.mat;
diff4=diff;
load entropydiff5.mat;
diff5=diff;
load entropydiff6.mat;
diff6=diff;
load entropydiff7.mat;
diff7=diff;
load entropydiff8.mat;
diff8=diff;
load entropydiff9.mat;
diff9=diff;
load entropydiff10.mat;
diff10=diff;
plot(1:1:length(diff2),diff2,'-*','LineWidth',2);
hold on;
plot(1:1:length(diff3),diff3,'LineWidth',2);
plot(1:1:length(diff4),diff4,'LineWidth',2);
plot(1:1:length(diff5),diff5,'LineWidth',2);
plot(1:1:length(diff6),diff6,'-+','LineWidth',2);
plot(1:1:length(diff7),diff7,'LineWidth',2);
plot(1:1:length(diff8),diff8,'LineWidth',2);
plot(1:1:length(diff9),diff9,'LineWidth',2);
plot(1:1:length(diff10),diff10,'-o','LineWidth',2);
legend('2x','3x','4x','5x','6x','7x','8x','9x','10x')
xlabel('Image sub-blocks')
ylabel('Entropy Differences')
title('Stego and Original Entropy Variation for Different Compression Rates')
grid on
grid minor
hold off


% figure
% subplot(2,2,1)
% % imshow(im)
% imshow(dicom,'DisplayRange',[])
% title('Original Cover')
% 
% subplot(2,2,2)
% % imshow(uint16(finalstego),'DisplayRange',[])
% montage(uint16(finalstego),'DisplayRange',[])
% % imshow(coverstego)
% caption=strcat('Stego Image at Compression: x',num2str(compressionRatio));
% title(caption)
% 
% subplot(2,2,3)
% imhist(im)
% title('Histogram of Cover Image')
% 
% subplot(2,2,4)
% % imhist(uint8(finalstego))
% imhist(uint16(finalstego))
% caption=strcat('Histogram of Stego Image at Compression: x',num2str(compressionRatio));
% title(caption)

% figure
% % montage(uint16(withTampterDetected),'DisplayRange',[])
% imshow(uint16(withTampterDetected))
% caption=strcat('Tampered Stego Image at Compression: x',num2str(compressionRatio));
% title(caption)
figure
imshowpair(uint16(finalstego),im,'montage')
title('Stego Vs Original Image')