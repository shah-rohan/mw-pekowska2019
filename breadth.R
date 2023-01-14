library(tidyverse)

ab8580 = read.table("ghists/AR16-4_AB-8580_tss50.ghist", header = T, stringsAsFactors = F)
emd07473 = read.table("ghists/AR16-4_EMD-07-473_tss50.ghist", header = T, stringsAsFactors = F)
tfpa5 = read.table("ghists/AR16-4_TF-PA5-40086_tss50.ghist", header = T, stringsAsFactors = F)

breadth_calc = function(ghistframe)
{
  breadthdf = data.frame(Gene = character(), upstream = numeric(), downstream = numeric(), bivalency = character())
  
  for (i in 1:length(ghistframe$Gene))
  {
    if (ghistframe[i,62]>=25)
    {
      upstream_pos = ((ghistframe[i,2:62] %>% detect_index(.f = ~ .x < 25, .dir = "backward", .default = 1))-60)*50
    }
    else {
      upstream_pos = ((ghistframe[i, 62:70] %>% detect_index(.f = ~ .x >= 25, .dir = "forward", .default = 1)) - 1)*50
    }
    
    if (ghistframe[i,70]>=25)
    {
      downstream_pos = (ghistframe[i, 70:122] %>% detect_index(.f = ~ .x < 25, .dir = "forward", .default = 122))*50+350
    }
    else {
      downstream_pos = ((ghistframe[i, 62:70] %>% detect_index(.f = ~ .x >= 25, .dir = "backward", .default = 1)) - 1)*50
    }
    
    breadthdf = rbind(breadthdf, data.frame(Gene = ghistframe$Gene[i], upstream = upstream_pos, downstream = downstream_pos))
    
    if (i %% 1000 == 0) {print(i)}
  }
  
  return(breadthdf)
}

breadth_ab8580 = breadth_calc(ab8580)
breadth_emd07473 = breadth_calc(emd07473)
breadth_tfpa5 = breadth_calc(tfpa5)

breadth_ab8580 = breadth_ab8580 %>% mutate(breadth = downstream - upstream)
breadth_emd07473 = breadth_emd07473 %>% mutate(breadth = downstream - upstream)
breadth_tfpa5 = breadth_tfpa5 %>% mutate(breadth = downstream - upstream)

bh4d = list()
bh4d$tf = (breadth_tfpa5 %>% filter(breadth >= 2000))$Gene
bh4d$ab = (breadth_ab8580 %>% filter(breadth >= 2000))$Gene
bh4d$emd = (breadth_emd07473 %>% filter(breadth >= 2000))$Gene

bh4d_ab_only = setdiff(bh4d$ab, bh4d$tf)
bh4d_emd_only = setdiff(bh4d$emd, bh4d$tf)
