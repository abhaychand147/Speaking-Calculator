#pragma once
#include "use.h"
#include "guidelines.h"
#include "train.h"

namespace Demo {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for menu
	/// </summary>
	public ref class menu : public System::Windows::Forms::Form
	{
	public:
		menu(void)
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
		~menu()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	protected: 
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::Button^  button2;
	private: System::Windows::Forms::Button^  button1;

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
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(menu::typeid));
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->BackColor = System::Drawing::SystemColors::ActiveCaptionText;
			this->label1->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 36, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label1->ForeColor = System::Drawing::Color::White;
			this->label1->Location = System::Drawing::Point(12, 37);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(587, 57);
			this->label1->TabIndex = 13;
			this->label1->Text = L"SPEAKING CALCULATOR";
			this->label1->Click += gcnew System::EventHandler(this, &menu::label1_Click);
			// 
			// button3
			// 
			this->button3->BackColor = System::Drawing::Color::Gray;
			this->button3->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 20.25F, System::Drawing::FontStyle::Bold));
			this->button3->Location = System::Drawing::Point(22, 357);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(539, 88);
			this->button3->TabIndex = 12;
			this->button3->Text = L"GUIDELINES";
			this->button3->UseVisualStyleBackColor = false;
			this->button3->Click += gcnew System::EventHandler(this, &menu::button3_Click);
			// 
			// button2
			// 
			this->button2->BackColor = System::Drawing::Color::Gray;
			this->button2->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 20.25F, System::Drawing::FontStyle::Bold));
			this->button2->Location = System::Drawing::Point(22, 247);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(539, 88);
			this->button2->TabIndex = 11;
			this->button2->Text = L"START CALCULATOR";
			this->button2->UseVisualStyleBackColor = false;
			this->button2->Click += gcnew System::EventHandler(this, &menu::button2_Click);
			// 
			// button1
			// 
			this->button1->BackColor = System::Drawing::Color::Gray;
			this->button1->Font = (gcnew System::Drawing::Font(L"Maiandra GD", 20.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->button1->Location = System::Drawing::Point(22, 137);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(539, 88);
			this->button1->TabIndex = 10;
			this->button1->Text = L"TRAIN";
			this->button1->UseVisualStyleBackColor = false;
			this->button1->Click += gcnew System::EventHandler(this, &menu::button1_Click);
			// 
			// menu
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
			this->ClientSize = System::Drawing::Size(1184, 761);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->button1);
			this->Name = L"menu";
			this->Text = L"menu";
			this->Load += gcnew System::EventHandler(this, &menu::menu_Load);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void menu_Load(System::Object^  sender, System::EventArgs^  e) {
			 }
private: System::Void label1_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
			train^ f2=gcnew train();
			 f2->ShowDialog(this);

		 }
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
			 use^ f2=gcnew use();
			 f2->ShowDialog(this);

		 }
private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) {
			 guidelines^ f3=gcnew guidelines();
			 f3->ShowDialog(this);
		 }
};
}
