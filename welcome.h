#pragma once
#include "menu.h"
namespace Demo {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for welcome
	/// </summary>
	public ref class welcome : public System::Windows::Forms::Form
	{
	public:
		welcome(void)
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
		~welcome()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label3;
	protected: 
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Label^  label1;

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
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(welcome::typeid));
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->SuspendLayout();
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->BackColor = System::Drawing::Color::Transparent;
			this->label3->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 24, System::Drawing::FontStyle::Bold));
			this->label3->ForeColor = System::Drawing::SystemColors::ButtonHighlight;
			this->label3->Location = System::Drawing::Point(350, 136);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(468, 78);
			this->label3->TabIndex = 11;
			this->label3->Text = L"THE SPEAKING CALCULATOR\r\n                CALCY";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->BackColor = System::Drawing::Color::Transparent;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft YaHei UI", 11.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label2->ForeColor = System::Drawing::SystemColors::ButtonHighlight;
			this->label2->Location = System::Drawing::Point(524, 107);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(104, 19);
			this->label2->TabIndex = 10;
			this->label2->Text = L"Welcome to \r\n";
			// 
			// button3
			// 
			this->button3->BackColor = System::Drawing::Color::Gray;
			this->button3->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 20.25F, System::Drawing::FontStyle::Bold));
			this->button3->ForeColor = System::Drawing::SystemColors::ButtonHighlight;
			this->button3->Location = System::Drawing::Point(343, 384);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(493, 76);
			this->button3->TabIndex = 9;
			this->button3->Text = L"GUIDELINES";
			this->button3->UseVisualStyleBackColor = false;
			this->button3->Click += gcnew System::EventHandler(this, &welcome::button3_Click);
			// 
			// button1
			// 
			this->button1->BackColor = System::Drawing::Color::Gray;
			this->button1->Font = (gcnew System::Drawing::Font(L"Rockwell", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->button1->ForeColor = System::Drawing::SystemColors::ButtonHighlight;
			this->button1->Location = System::Drawing::Point(460, 588);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(283, 65);
			this->button1->TabIndex = 8;
			this->button1->Text = L"Continue to Calcy ->";
			this->button1->UseVisualStyleBackColor = false;
			this->button1->Click += gcnew System::EventHandler(this, &welcome::button1_Click);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->BackColor = System::Drawing::Color::Transparent;
			this->label1->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 20.25F, System::Drawing::FontStyle::Bold));
			this->label1->ForeColor = System::Drawing::SystemColors::ButtonHighlight;
			this->label1->Location = System::Drawing::Point(312, 255);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(567, 96);
			this->label1->TabIndex = 7;
			this->label1->Text = L"KINDLY GO THROUGH THE GUIDELINES \r\n          TO USE THE CALCULATOR\r\n             " 
				L"BEFORE YOU PROCEED";
			// 
			// welcome
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
			this->ClientSize = System::Drawing::Size(1184, 761);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->label1);
			this->Name = L"welcome";
			this->Text = L"welcome";
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion


private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) {

			 guidelines^ f3=gcnew guidelines();
			 f3->ShowDialog(this);
		 }
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
			 menu^ f3=gcnew menu();
			 this->Hide();

			 f3->ShowDialog(this);
		 }
};
}
