#ifndef MAINLOOP_HELPER_H_
#define MAINLOOP_HELPER_H_

#include <csl.h>
#include <csl_mcbsp.h>
#include <csl_edma.h>
#include "enkoder.h"

//#include "ADC8361.h"		/*ADC utility data*/
#include "DAC7716.h"		/*DAC utility data*/

#define RESET 0x300000
#define ACTIVATE 0x700000
#define R_HIGH 0x200000
#define R_LOW 0x600000

// CUDA POTRZEBNE DO DZIA�ANIA EMIF
#define ADDRESS 0xA0000000	// adres zapisu z EMIF CE2 (Daughter Card) - tu wystawiamy wyjscia
#define GBLCTL 0x01800000 // nie wiem po co to....
#define CE2 0x01800010	//

/*Constants*/
static unsigned short int const ML_BUTTON1 = 0x1;
static unsigned short int const ML_BUTTON2 = 0x2;
static unsigned short int const ML_BUTTON3 = 0x4;
static unsigned short int const ML_BUTTON4 = 0x8;

#define TRUE 1
#define FALSE 0

/*Variables*/
/*Flags*/
static volatile short Enc_ReceiveToGo = TRUE;
static volatile short DAC_TransmitToGo = TRUE;

/*Prototypes*/
void ML_ExitProcedure(EDMA_Handle* _hEDMA_ADC, int _ADC_IRQ, EDMA_Handle* _hEDMA_ADC_Ping, EDMA_Handle* _hEDMA_ADC_Pong, EDMA_Handle* _hEDMA_Dac, int _DAC_IRQ, EDMA_Handle* _hEDMA_DAC_Ping,  EDMA_Handle* _hEDMA_DAC_Pong,  MCBSP_Handle _hMCBSP_ADC, MCBSP_Handle _hMCBSP_DAC, int _ShowMsg);
unsigned short int ML_CheckSwitch();
//ADC_Measurement ML_WaitForFreshInput(volatile ADC_Measurement* _ADCInput, volatile short* _ReceiveToGo);
Enc_Measurement Enc_WaitForFreshInput(volatile Enc_Measurement* _EncInput, volatile short*_ReceiveToGo);
DAC_Values ML_PrepareFreshOutput( volatile short* _TransmitToGo);
DAC_Values Enc_PrepareFreshOutput(Enc_Measurement _EncInput, volatile short* _TransmitToGo);
void Enc_SendOrder(int order);


#endif /* MAINLOOP_HELPER_H_ */
