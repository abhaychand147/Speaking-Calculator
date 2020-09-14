#pragma once

#include "iostream"
#include "train_test_module.h"
#include "speakTheResult.h"

using namespace std;
namespace Demo {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for use
	/// </summary>
	public ref class use : public System::Windows::Forms::Form
	{
	public:
		use(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~use()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Button^  button2;
	protected: 
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label6;

	protected: 





	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(use::typeid));
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->SuspendLayout();
			// 
			// button2
			// 
			this->button2->BackColor = System::Drawing::Color::Black;
			this->button2->Font = (gcnew System::Drawing::Font(L"Rockwell", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->button2->ForeColor = System::Drawing::Color::Olive;
			this->button2->Location = System::Drawing::Point(18, 317);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(305, 65);
			this->button2->TabIndex = 11;
			this->button2->Text = L"<- BACK";
			this->button2->UseVisualStyleBackColor = false;
			this->button2->Click += gcnew System::EventHandler(this, &use::button2_Click);
			// 
			// button1
			// 
			this->button1->BackColor = System::Drawing::Color::Black;
			this->button1->Font = (gcnew System::Drawing::Font(L"Rockwell", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->button1->ForeColor = System::Drawing::Color::Olive;
			this->button1->Location = System::Drawing::Point(18, 246);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(305, 65);
			this->button1->TabIndex = 10;
			this->button1->Text = L"GO CALCY ->";
			this->button1->UseVisualStyleBackColor = false;
			this->button1->Click += gcnew System::EventHandler(this, &use::button1_Click);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->BackColor = System::Drawing::Color::Transparent;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei", 18, System::Drawing::FontStyle::Italic, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label1->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->label1->Location = System::Drawing::Point(12, 56);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(1068, 217);
			this->label1->TabIndex = 9;
			this->label1->Text = resources->GetString(L"label1.Text");
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->BackColor = System::Drawing::Color::Transparent;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 36, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label5->ForeColor = System::Drawing::Color::FloralWhite;
			this->label5->Location = System::Drawing::Point(91, 524);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(222, 64);
			this->label5->TabIndex = 16;
			this->label5->Text = L"EQUALS";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->BackColor = System::Drawing::Color::Transparent;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 36, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label4->ForeColor = System::Drawing::Color::Cyan;
			this->label4->Location = System::Drawing::Point(101, 597);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(179, 64);
			this->label4->TabIndex = 15;
			this->label4->Text = L"Result";
			this->label4->Click += gcnew System::EventHandler(this, &use::label4_Click);
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->BackColor = System::Drawing::Color::Transparent;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 20.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label3->ForeColor = System::Drawing::Color::DarkOliveGreen;
			this->label3->Location = System::Drawing::Point(165, 468);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(134, 36);
			this->label3->TabIndex = 14;
			this->label3->Text = L"operator";
			this->label3->Click += gcnew System::EventHandler(this, &use::label3_Click);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->BackColor = System::Drawing::Color::Transparent;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 36, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label2->ForeColor = System::Drawing::Color::FloralWhite;
			this->label2->Location = System::Drawing::Point(320, 449);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(129, 64);
			this->label2->TabIndex = 13;
			this->label2->Text = L"Op2";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->BackColor = System::Drawing::Color::Transparent;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 36, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label6->ForeColor = System::Drawing::Color::FloralWhite;
			this->label6->Location = System::Drawing::Point(12, 449);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(129, 64);
			this->label6->TabIndex = 12;
			this->label6->Text = L"Op1";
			this->label6->Click += gcnew System::EventHandler(this, &use::label6_Click);
			// 
			// use
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
			this->ClientSize = System::Drawing::Size(1184, 761);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->label1);
			this->Name = L"use";
			this->Text = L"use";
			this->Load += gcnew System::EventHandler(this, &use::use_Load);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void label1_Click(System::Object^  sender, System::EventArgs^  e) {

			 }
private: System::Void use_Load(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {

	//cout<<"Insert Testing module";	
			// step5();
			 int op1 = 2, op2 = 3;
			 char opt = '/';
			int res = 0;
			
			
				step5(&op1,&opt,&op2);


				
				
				//string ^ in;

				String ^ in1=" ";
				String ^ in2=" ";
				String ^ in3=" ";
				String ^ in4=" ";

				//^in2,^in3,^in4;
					in1=System::Convert::ToString(op1);
				//String ^ in2;
					in2=System::Convert::ToString(op2);
				//String ^ in3;
					//in3=System::Convert::ToString(opt);
				//String ^ in4;
					//in4=System::Convert::ToString(res);

				switch(opt) {
				case '/': // by
					
						res = op1 / op2;
						//in="/";
						label3->Text="/";
						 label6->Text=in1;
						 //label3->Text=in;
						 label2->Text=in2;
						 in4=System::Convert::ToString(res);
						 label4->Text=in4;
						speakResult(res);
						break;
				case '*' : // into
						res = op1 * op2;
						//in="*";
						label3->Text="*";
						label6->Text=in1;
						 //label3->Text=in;
						 label2->Text=in2;
						 in4=System::Convert::ToString(res);
						 label4->Text=in4;
						speakResult(res);
						break;
				case '-' : // minus
						res = op1 - op2;
						//in="-";
						label3->Text="-";
						 label6->Text=in1;
						 //label3->Text=in;
						 label2->Text=in2;
						 in4=System::Convert::ToString(res);
						 label4->Text=in4;
						speakResult(res);
						break;
				case '+' : // plus
						res = op1 + op2;
						//in="+";
						label3->Text="+";
						 label6->Text=in1;
						 //label3->Text=in;
						 label2->Text=in2;
						 in4=System::Convert::ToString(res);
						 label4->Text=in4;
						speakResult(res);
						break;
				default :
						break;
				}





	//--------------------------------------------------------------
	// Displaying the result


		 }
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
			 this->Hide();
		 }
private: System::Void label3_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label4_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void label6_Click(System::Object^  sender, System::EventArgs^  e) {
 
		 }
};
}
