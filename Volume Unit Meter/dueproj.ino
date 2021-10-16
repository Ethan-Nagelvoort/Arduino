/*
 * FILENAME : VUMeterFirmware.ino
 * DATE     : June 7th, 2021
 * 
 * AUTHOR 1 : Robert Esposito
 * AUTHOR 2 : Ethan Nagelvoort
 * COURSE   : EE/CompE 496B - Senior Design B
 * TEAM     : Team 10 - VU Meter
 * 
 * DESCRIPTION :
 *     This program will read from the Arduino Due's ADC, convert those values to
 *     a decible scale, and then big bang the results to shift registers.
 * 
 * LAST MODIFIED : 10/10/2021
 * 
 */

volatile int clockCount = 17;
volatile int dataCount = 7;
int a[8];
uint8_t n = 0;
uint8_t data[8] = {0,0,0,0,0,0,0,0};
uint8_t temp[8];


void setup(){
  /**************** ADC SETTINGS ****************/
  ADC->ADC_MR |= ADC_MR_FREERUN | ADC_MR_LOWRES;
  ADC->ADC_CR = ADC_CR_START;
  ADC->ADC_CHER = ADC_CHER_CH7 | ADC_CHER_CH6 | ADC_CHER_CH5 | ADC_CHER_CH4 | ADC_CHER_CH3 | ADC_CHER_CH2 | ADC_CHER_CH1 | ADC_CHER_CH0;

  /**************** OUTPUT PIN SETTINGS ****************/
  //PA15 (Pin 24) = CLOCK.
  PIOA->PIO_OER |= PIO_OER_P15;
  PIOA->PIO_OWER |= PIO_OWER_P15;
  
  //PC22 (Pin 8) = LATCH.
  PIOC->PIO_OER |= PIO_OER_P22;
  PIOC->PIO_OWER |= PIO_OWER_P22;

  //PD7 (Pin 11) = Shift Register 0.
  PIOD->PIO_OER |= PIO_OER_P7;
  PIOD->PIO_OWER |= PIO_OWER_P7;

  //PD8 (Pin 12) = Shift Register 1.
  PIOD->PIO_OER |= PIO_OER_P8;
  PIOD->PIO_OWER |= PIO_OWER_P8;

  //PB25 (Pin 2) = Shift Register 2.
  PIOB->PIO_OER |= PIO_OER_P25;
  PIOB->PIO_OWER |= PIO_OWER_P25;
  
  //PC28 (Pin 3) = Shift Register 3.
  PIOC->PIO_OER |= PIO_OER_P28;
  PIOC->PIO_OWER |= PIO_OWER_P28;
  
  //PC26 (Pin 4) = Shift Register 4.
  PIOC->PIO_OER |= PIO_OER_P26;
  PIOC->PIO_OWER |= PIO_OWER_P26;
  
  //PC25 (Pin 5) = Shift Register 5.
  PIOC->PIO_OER |= PIO_OER_P25;
  PIOC->PIO_OWER |= PIO_OWER_P25;
  
  //PC24 (Pin 6) = Shift Register 6.
  PIOC->PIO_OER |= PIO_OER_P24;
  PIOC->PIO_OWER |= PIO_OWER_P24;
  
  //PC23 (Pin 7) = Shift Register 7.
  PIOC->PIO_OER |= PIO_OER_P23;
  PIOC->PIO_OWER |= PIO_OWER_P23;

  /**************** TIMER AND INTERRUPT SETTINGS ****************/
  PMC->PMC_PCR = PMC_PCR_EN | PMC_PCR_CMD | (ID_TC0 & 0x7F);
  TC0->TC_CHANNEL[0].TC_CMR = TC_CMR_WAVE | TC_CMR_WAVSEL_UP_RC | TC_CMR_TCCLKS_TIMER_CLOCK1;
  TC0->TC_CHANNEL[0].TC_RC = 233;
  TC0->TC_CHANNEL[0].TC_IER = TC_IER_CPCS;
  NVIC_EnableIRQ(TC0_IRQn);
  TC0->TC_CHANNEL[0].TC_CCR = TC_CCR_CLKEN | TC_CCR_SWTRG;

  /**************** CLOCK START ****************/
  PIOA->PIO_ODSR |= PIO_ODSR_P15;
}

void loop(){
  if (clockCount < 0){
    clockCount = 17;
  }
  
  if (dataCount < 0){
    dataCount = 7;
    a[0]=ADC->ADC_CDR[7];   // a[0] (analog pin 0) = sub bass filter.
    a[1]=ADC->ADC_CDR[6];   // a[1] (analog pin 1) = bass filter.
    a[2]=ADC->ADC_CDR[5];   // a[2] (analog pin 2) = lower midrange filter.
    a[3]=ADC->ADC_CDR[4];   // a[3] (analog pin 3) = midrange filter.
    a[4]=ADC->ADC_CDR[3];   // a[4] (analog pin 4) = upper midrange filter.
    a[5]=ADC->ADC_CDR[2];   // a[5] (analog pin 5) = presence filter.
    a[6]=ADC->ADC_CDR[1];   // a[6] (analog pin 6) = brilliance filter.
    a[7]=ADC->ADC_CDR[0];   // a[7] (analog pin 7) = entire frequency filter.
    
    for(int i = 0; i < 8; i++){
      if(a[i] < 77) {
        n=0;                             //0 LEDs.
      }
      else if(a[i] < 109){
        n=1;                             //1 LEDs.
      }
      else if(a[i] < 154){
        n=2;                             //2 LEDs.
      }
      else if(a[i] < 217){
        n=3;                             //3 LEDs.
      }
      else if(a[i] < 306){
        n=4;                             //4 LEDs.
      }
      else if(a[i] < 432){
        n=5;                             //5 LEDs.
      }
      else if(a[i] < 610){
        n=6;                             //6 LEDs.
      }
      else if(a[i] < 861){
        n=7;                             //7 LEDs.
      }
      else{
        n=8;                             //8 LEDs.
      }
      data[i] = 0xFF >> 8 - n;
    }
  }
}

void TC0_Handler(){

  TC0->TC_CHANNEL[0].TC_SR;
  PIOA->PIO_ODSR ^= PIO_ODSR_P15;
  
  if (clockCount > 2){
    PIOC->PIO_ODSR &= ~PIO_ODSR_P22;
    
    if ((clockCount & 0x01) != 0){
      for (int i = 0; i < 8; i++){
        temp[i] = data[i] & (0x01 << dataCount);
        temp[i] >>= dataCount;
      }
      dataCount--;

      //Shift Register 0.
      if (temp[0] == 0x01){
        PIOD->PIO_ODSR |= PIO_ODSR_P7;
      }
      else{
        PIOD->PIO_ODSR &= ~PIO_ODSR_P7;
      }
      
      //Shift Register 1.
      if (temp[1] == 0x01){
        PIOD->PIO_ODSR |= PIO_ODSR_P8;
      }
      else{
        PIOD->PIO_ODSR &= ~PIO_ODSR_P8;
      }

      //Shift Register 2.
      if (temp[2] == 0x01){
        PIOB->PIO_ODSR |= PIO_ODSR_P25;
      }
      else{
        PIOB->PIO_ODSR &= ~PIO_ODSR_P25;
      }

      //Shift Register 3.
      if (temp[3] == 0x01){
        PIOC->PIO_ODSR |= PIO_ODSR_P28;
      }
      else{
        PIOC->PIO_ODSR &= ~PIO_ODSR_P28;
      }

      //Shift Register 4.
      if (temp[4] == 0x01){
        PIOC->PIO_ODSR |= PIO_ODSR_P26;
      }
      else{
        PIOC->PIO_ODSR &= ~PIO_ODSR_P26;
      }   

      //Shift Register 5.
      if (temp[5] == 0x01){
        PIOC->PIO_ODSR |= PIO_ODSR_P25;
      }
      else{
        PIOC->PIO_ODSR &= ~PIO_ODSR_P25;
      } 

      //Shift Register 6.
      if (temp[6] == 0x01){
        PIOC->PIO_ODSR |= PIO_ODSR_P24;
      }
      else{
        PIOC->PIO_ODSR &= ~PIO_ODSR_P24;
      } 

      //Shift Register 7.
      if (temp[7] == 0x01){
        PIOC->PIO_ODSR |= PIO_ODSR_P23;
      }
      else{
        PIOC->PIO_ODSR &= ~PIO_ODSR_P23;
      }
    }
  }
  else{
    PIOC->PIO_ODSR |= PIO_ODSR_P22;
  }
  clockCount--;
}
